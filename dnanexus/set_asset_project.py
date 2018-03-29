#!/usr/bin/env python

from __future__ import print_function
import sys
import dxpy
import json
from collections import OrderedDict

try:
    project = dxpy.DXProject(sys.argv[1])
except IndexError:
    print("Must supply a DNAnexus project name or ID", file=sys.stderr)
    sys.exit(1)
except dxpy.DXError:
    projects = dxpy.find_projects(sys.argv[1])
    try:
        project_id = next(projects).get('id')
    except StopIteration:
        print("Could not resolve %s to a DNAnexus project." % (sys.argv[1]), file=sys.stderr)
        sys.exit(1)
    else:
        project = dxpy.DXProject(project_id)

applets = [
    'bam2tagalign',
    'encode_map',
    'filter_qc',
    'xcor',
    'xcor_only',
    'spp',
    'pool',
    'pseudoreplicator',
    'encode_spp',
    'encode_macs2',
    'macs2',
    'idr2',
    'encode_idr',
    'overlap_peaks',
    'input_shield',
    'accession_analysis',
    'shell',
    'scrub'
]

for applet in applets:
    fh = open("%s/dxapp.json" % (applet), 'r')
    try:
        dxapp_json = json.load(fh, object_pairs_hook=OrderedDict)
    except ValueError:
        print("Cannot interpret JSON in %s/dxapp.json" % (applet))
        sys.exit(1)
    asset_depends = dxapp_json['runSpec']['assetDepends']
    for asset in asset_depends:
        asset.update(project=project.get_id())
    print(json.dumps(dxapp_json, indent=4, separators=(',', ': ')))
