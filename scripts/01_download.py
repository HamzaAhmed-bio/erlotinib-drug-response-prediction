"""Download frozen public Sanger inputs and record hashes; no credentials needed."""
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import urllib.request
import io
import struct
import zipfile
import zlib
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get('DRUG_DATA_DIR', ROOT / 'data'))
FILES = {
    'GDSC2_fitted_dose_response_27Oct23.xlsx': 'https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx',
    'model_list_20241120.csv': 'https://cog.sanger.ac.uk/cmp/download/model_list_20241120.csv',
    'GDSC_Fitted_Data_Description.pdf': 'https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC_Fitted_Data_Description.pdf',
}
ARCHIVE = 'https://cog.sanger.ac.uk/cmp/download/rnaseq_merged_Sanger_Broad_20210329.zip'
MEMBER = 'rnaseq_fpkm_20210329.csv'

class RemoteZip(io.RawIOBase):
    def __init__(self, url):
        self.url = url
        with urllib.request.urlopen(urllib.request.Request(url, method='HEAD'), timeout=60) as r:
            self.size = int(r.headers['Content-Length'])
            self.etag = r.headers.get('ETag')
        self.pos = 0
    def seek(self, offset, whence=0):
        self.pos = offset if whence == 0 else self.pos+offset if whence == 1 else self.size+offset
        return self.pos
    def tell(self): return self.pos
    def read(self, size=-1):
        if size < 0: size = self.size-self.pos
        if size == 0: return b''
        req = urllib.request.Request(self.url, headers={'Range':f'bytes={self.pos}-{self.pos+size-1}'})
        with urllib.request.urlopen(req, timeout=120) as r:
            assert r.status == 206, 'Server did not honor byte range'
            result = r.read()
        assert len(result) == size
        self.pos += size
        return result

def download_expression():
    path = DATA / MEMBER
    remote = RemoteZip(ARCHIVE)
    with zipfile.ZipFile(remote) as z:
        info = z.getinfo(MEMBER)
    if not path.exists():
        remote.seek(info.header_offset)
        header = remote.read(30)
        name_len, extra_len = struct.unpack('<HH', header[26:30])
        remote.seek(info.header_offset+30+name_len+extra_len)
        print('Downloading FPKM member only:', info.compress_size, 'compressed bytes', flush=True)
        packed = remote.read(info.compress_size)
        assert info.compress_type == zipfile.ZIP_DEFLATED
        unpacked = zlib.decompress(packed, -15)
        assert len(unpacked) == info.file_size
        assert zlib.crc32(unpacked) == info.CRC
        path.write_bytes(unpacked)
    with path.open('rb') as f:
        digest = hashlib.file_digest(f, 'sha256').hexdigest()
    crc = 0
    with path.open('rb') as f:
        while block := f.read(1024*1024): crc = zlib.crc32(block, crc)
    assert crc == info.CRC and path.stat().st_size == info.file_size
    print(MEMBER, path.stat().st_size, digest, flush=True)
    return {'file':MEMBER, 'url':ARCHIVE, 'archive_member':MEMBER, 'archive_etag':remote.etag,
            'archive_bytes':remote.size, 'zip_crc32_verified':True, 'bytes':path.stat().st_size, 'sha256':digest}

def download(item):
    name, url = item
    path = DATA / name
    if not path.exists():
        request = urllib.request.Request(url, headers={'User-Agent':'ResearchPortfolio/1.0'})
        with urllib.request.urlopen(request, timeout=120) as response, path.with_suffix(path.suffix+'.part').open('wb') as out:
            print(name, 'downloading', response.headers.get('Content-Length','unknown'), flush=True)
            while block := response.read(1024*1024):
                out.write(block)
        path.with_suffix(path.suffix+'.part').replace(path)
    digest = hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()
    print(name, path.stat().st_size, digest, flush=True)
    return {'file': name, 'url':url, 'bytes':path.stat().st_size, 'sha256':digest}

if __name__ == '__main__':
    DATA.mkdir(parents=True, exist_ok=True)
    (ROOT/'metadata').mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(download, FILES.items()))
    results.append(download_expression())
    receipt = {'retrieved_utc':datetime.now(timezone.utc).isoformat(), 'files':results}
    receipt_path = ROOT/'metadata/source_receipt.json'
    if receipt_path.exists():
        old = json.loads(receipt_path.read_text())
        assert {x['file']:x['sha256'] for x in old['files']} == {x['file']:x['sha256'] for x in results}, 'Upstream data changed'
    else:
        receipt_path.write_text(json.dumps(receipt, indent=2))
