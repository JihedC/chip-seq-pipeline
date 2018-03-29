#!/usr/bin/env python

from __future__ import print_function
import sys
import dxpy
import json
from collections import OrderedDict
import shutil

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
    infile_path = "%s/dxapp.json" % (applet)
    bufile_path = "%s/dxapp-bu.json" % (applet)
    outfile_path = "%s/dxapp.json" % (applet)
    infh = open(infile_path, 'r')
    try:
        dxapp_json = json.load(infh, object_pairs_hook=OrderedDict)
    except ValueError:
        print("Cannot interpret JSON in %s/dxapp.json" % (applet), file=sys.stderr)
        sys.exit(1)
    asset_depends = dxapp_json['runSpec']['assetDepends']
    for asset in asset_depends:
        asset.update(project=project.get_id())
    # test to make sure serialization works before destroying the dxapp.json file
    try:
        s = json.dumps(dxapp_json, indent=4, separators=(',', ': '))
    except:
        print("Cannot serialize new JSON", file=sys.stderr)
        sys.exit(1)
    infh.close()
    shutil.copyfile(infile_path, bufile_path)
    out_fh = open(outfile_path, 'w+')
    json.dump(dxapp_json, out_fh, indent=4, separators=(',', ': '))
