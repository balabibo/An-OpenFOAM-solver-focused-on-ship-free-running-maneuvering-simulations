# An OpenFOAM solver focused on free-running maneuvering simulations

Summary

These modifications are based on the "overInterDyMFoam" solver in OpenFOAM v2312.
There are several modules which were modified or added, and its general function is realizing the maneuvering motions of ship.
A 35 degree turning tutorial has been updated.

Functions
  1. Modified the storage structure of rigidBodyState dictionary, which can make change of number of Degree of freedom (Dof) when you modify the dynamicMeshDict in /constant possible.
  2. add two solidBody motions "driven3DofMotion/driven2DofMotion" in /src/dynamicMesh/motionSolvers/displacement/solidBody/solidBodyMotionFunctions. It was modified from drivenLinearMotion, and it can realize the following rotation of background mesh region, which is useful for maneuvering motions like turning or zigzag.
  3. add two momentum source methods "oumSource/HOSource" in /src/fvOptions/sources/derived. They are based on blade element momentum theory/open water curve seperately, and you can find details from the reference, it can be used to replace the real propeller as a propulsion device.
  4. add 4 different maneuvering motions, self-propulsion, turning, zigzag and coursekeeping, for ship. First, it utilize the PID controller to adjust the revolution speed of discretized propeller or momentum source in "sailing" mode, and the PID contorller is also applied to control the rudder motion in "coursekeeping" mode. The rudder controller is used to realize "turning" and "zigzag" maneuvering motions.

P.S.

This serial will be updated continuously in the future. Those modules may have some bugs due to the neglect of author, and if you have any question or suggestion, please feel free to use the "Issues" button.
