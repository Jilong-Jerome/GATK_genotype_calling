#!/usr/bin/env python3
import sys, os, glob

sys.path.insert(0, os.path.realpath(
    os.path.join(os.path.dirname(__file__), 'workflow_sources')
))

from gwf import Workflow
from workflow_sources import gatk_calling_workflow

gwf = Workflow()

# Every *.config.yaml under configurations/ is loaded and registered as a set of
# targets. The shipped config.template.yaml is intentionally NOT matched by this
# glob, so it is never executed; copy it to configurations/<name>.config.yaml first.
configs = glob.glob(os.path.join(os.path.dirname(__file__), 'configurations', '*.config.yaml'))
for config in configs:
    gwf = gatk_calling_workflow(config_file=config, gwf=gwf)
