/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "HelloWorldAlg.h"
#include "GaudiKernel/MsgStream.h"

DECLARE_COMPONENT(HelloWorldAlg)

HelloWorldAlg::HelloWorldAlg(const std::string& aName, ISvcLocator* aSvcLoc) : Gaudi::Algorithm(aName, aSvcLoc) {}

HelloWorldAlg::~HelloWorldAlg() {}

StatusCode HelloWorldAlg::initialize() {
  if (Gaudi::Algorithm::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

StatusCode HelloWorldAlg::execute(const EventContext& ctx) const {
  info() << endmsg;
  info() << endmsg;
  info() << theMessage << endmsg;
  info() << endmsg;
  info() << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode HelloWorldAlg::finalize() { return Gaudi::Algorithm::finalize(); }
