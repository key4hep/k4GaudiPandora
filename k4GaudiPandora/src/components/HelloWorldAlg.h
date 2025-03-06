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

// Note that in most cases it's better to use the functional
// approach than Gaudi::Algorithm

// GAUDI
#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"

class HelloWorldAlg : public Gaudi::Algorithm {
public:
  HelloWorldAlg(const std::string&, ISvcLocator*);
  virtual ~HelloWorldAlg();
  /**  Initialize.
   *   @return status code
   */
  StatusCode initialize() override;
  /**  Execute.
   *   @return status code
   */
  StatusCode execute(const EventContext&) const override;
  /**  Finalize.
   *   @return status code
   */
  StatusCode finalize() override;

private:
  // member variable
  Gaudi::Property<std::string> theMessage{this, "PerEventPrintMessage", "Hello ",
                                          "The message to printed for each Event"};
};
