#!/usr/bin/env python
#
# Copyright (c) 2020-2024 Key4hep-Project.
#
# This file is part of Key4hep.
# See https://key4hep.github.io/key4hep-doc/ for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

# A simple script to compare the hits from the Gaudi and Marlin output
import argparse
from podio.root_io import Reader

parser = argparse.ArgumentParser(description="Compare hits from Gaudi and Marlin")
parser.add_argument(
    "--gaudi-file", default="output_calo_digi.root", help="Gaudi output file"
)
parser.add_argument(
    "--marlin-file", default="output_REC.edm4hep.root", help="Marlin output file"
)

parser.add_argument(
    "--gaudi-hits",
    default="GaudiRelationCaloHitECALBarrel",
    help="Gaudi hits collection",
)
parser.add_argument(
    "--marlin-hits", default="RelationCaloHitECALBarrel", help="Marlin hits collection"
)

args = parser.parse_args()

reader_gaudi = Reader(args.gaudi_file)
reader_marlin = Reader(args.marlin_file)

events_gaudi = reader_gaudi.get("events")
events_marlin = reader_marlin.get("events")

for i, frame_gaudi in enumerate(events_gaudi):
    frame_marlin = events_marlin[i]
    rels_gaudi = frame_gaudi.get(args.gaudi_hits)
    rels_marlin = frame_marlin.get(args.marlin_hits)
    print(f"Checking event {i} with {len(rels_gaudi)} hits")
    assert len(rels_gaudi) == len(
        rels_marlin
    ), f"Number of relations differ: {len(rels_gaudi)} vs {len(rels_marlin)}"
    for j, (rel_gaudi, rel_marlin) in enumerate(zip(rels_gaudi, rels_marlin)):
        print(f"Checking rel {j}")
        assert (
            rel_gaudi.getFrom() == rel_marlin.getFrom()
        ), f"From differs for relation {j}: {rel_gaudi.getFrom()} vs {rel_marlin.getFrom()}"
