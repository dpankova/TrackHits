#
# $Id: CMakeLists.txt 141949 2016-04-27 16:53 kjero $
#

i3_project(StartingTrackVeto
  PYTHON_DIR python)

i3_add_library(StartingTrackVeto
  private/StartingTrackVeto/StartingTrackVeto.cxx    
  private/StartingTrackVeto/StartingTrackVetoUtils.cxx    
  USE_PROJECTS icetray dataclasses recclasses phys-services
               photonics-service serialization
  USE_TOOLS boost python
)
