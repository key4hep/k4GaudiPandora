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
    "--gaudi-collections",
    default=[
        "GaudiECALBarrel",
        "GaudiECALEndcap",
        "GaudiHCALBarrel",
        "GaudiHCALEndcap",
        "GaudiHCALRing",
    ], nargs="+",
    help="Gaudi hits collection",
)
parser.add_argument(
    "--marlin-collections",
    default=[
        "ECALBarrel",
        "ECALEndcap",
        "HCALBarrel",
        "HCALEndcap",
        "HCALOther",
    ], nargs="+", help="Marlin hits collection"
)

args = parser.parse_args()

reader_gaudi = Reader(args.gaudi_file)
reader_marlin = Reader(args.marlin_file)

events_gaudi = reader_gaudi.get("events")
events_marlin = reader_marlin.get("events")

for i, frame_gaudi in enumerate(events_gaudi):
    frame_marlin = events_marlin[i]
    for collection_marlin, collection_gaudi in zip(args.marlin_collections, args.gaudi_collections):
        print(f'Checking collection "{collection_marlin}" (Marlin) and "{collection_gaudi}" (Gaudi)')
        hits_gaudi = frame_gaudi.get(collection_gaudi)
        hits_marlin = frame_marlin.get(collection_marlin)
        print(f"Checking event {i} with {len(hits_gaudi)} hits")
        assert len(hits_gaudi) == len(
            hits_marlin
        ), f"Number of hits differ: {len(hits_gaudi)} vs {len(hits_marlin)}"
        for j, (hit_gaudi, hit_marlin) in enumerate(zip(hits_gaudi, hits_marlin)):
            print(f"Checking hit {j}")
            assert (
                hit_gaudi.getCellID() == hit_marlin.getCellID()
            ), f"CellID differ for hit {j}: {hit_gaudi.getCellID()} vs {hit_marlin.getCellID()}"
            assert (
                hit_gaudi.getEnergy() == hit_marlin.getEnergy()
            ), f"Energy differ for hit {j}: {hit_gaudi.getEnergy()} vs {hit_marlin.getEnergy()}"
            assert (
                hit_gaudi.getEnergyError() == hit_marlin.getEnergyError()
            ), f"EnergyError differ for hit {j}: {hit_gaudi.getEnergyError()} vs {hit_marlin.getEnergyError()}"
            assert (
                hit_gaudi.getTime() == hit_marlin.getTime()
            ), f"Time differ for hit {j}: {hit_gaudi.getTime()} vs {hit_marlin.getTime()}"
            assert (
                hit_gaudi.getPosition() == hit_marlin.getPosition()
            ), f"Position differ for hit {j}: {hit_gaudi.getPosition()} vs {hit_marlin.getPosition()}"
            assert (
                hit_gaudi.getType() == hit_marlin.getType()
            ), f"Type differ for hit {j}: {hit_gaudi.getType()} vs {hit_marlin.getType()}"
