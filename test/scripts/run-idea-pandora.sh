#!/bin/bash
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

# Run the IDEA Pandora job for the CI.
#
# PandoraSettingsIdea.xml refers to its ONNX models by bare file name, i.e. relative to
# the working directory, so link the copies that ExternalData fetched into the build tree
# next to the job.  The remaining arguments are passed on to k4run.
#
# Usage: run-idea-pandora.sh <model.onnx> <model.onnx> <options.py> [k4run args...]
set -e

ln -sf "$1" "$(basename "$1")"
ln -sf "$2" "$(basename "$2")"
shift 2

exec k4run "$@"
