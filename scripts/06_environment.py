"""Record installed versions for the dependency closure of this project."""
import importlib.metadata as metadata
import json
import platform
import sys
from pathlib import Path
from packaging.requirements import Requirement
from packaging.utils import canonicalize_name

root=Path(__file__).resolve().parents[1]
queue=[Requirement(x).name for x in (root/'requirements.txt').read_text().splitlines()]
found={}
while queue:
    name=queue.pop()
    key=canonicalize_name(name)
    if key in found: continue
    distribution=metadata.distribution(name)
    found[key]=distribution.version
    for spec in distribution.requires or []:
        requirement=Requirement(spec)
        if requirement.marker is None or requirement.marker.evaluate({'extra':''}):
            queue.append(requirement.name)
(root/'requirements-lock.txt').write_text('\n'.join(k+'=='+v for k,v in sorted(found.items()))+'\n')
(root/'metadata/environment.json').write_text(json.dumps({'python':sys.version,'platform':platform.platform(),'packages':found},indent=2))
print(len(found),'pinned packages')
