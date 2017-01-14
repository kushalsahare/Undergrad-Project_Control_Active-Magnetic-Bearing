# Performance Analysis Design and Control of Rotors Mounted with Active Magnetic Bearings

This is my Undergraduate thesis work. Guided by Prof Rajiv Tiwari and Prof Atanu Banerjee, department of Mechanical Engineering at IIT Guwahati, India.

Abstract

Advanced rotor systems today consist of a lightweight rotor supported by radial active magnetic bearings
(AMB). These systems are widely used in flywheel applications and in other fields where high rotational
speeds are essential. Unbalance response is a common vibration problem associated with rotating
machinery. For several years, researchers have demonstrated that this vibration could be greatly alleviated
for machines using active magnetic bearings through active control.
Monitoring of response of rotors and handling the unbalance in the rotors is a high priority task in the
field of rotor dynamics. Various modeling techniques have been used for this purpose some of which
include modeling of Euler Bernoulli and Timoshenko beams with various subsequent modeling
techniques like Modal analysis, energy methods(using Galerkin analysis), Simulink etc.
In this paper work has been done on modeling of rotor as Euler Bernoulli beams with AMBs. FEM (Finite
Element Method) has been used to model the rotor-bearing system which includes variations such as shaft
as massless and with mass, rigid and flexible, single disk or two disk rotor system, single AMB or multi
AMB system, under the influence of constant forces and also impact forces at selected locations. The
complexity of the model has been gradually increased from 3 element model to 9 element model.
Further on the response magnitude and phase for all the above cases have been plotted with spin speeds
and critical points noted with and without the action of AMBs. Then the unbalance problem in rotors has
been tackled by again gradually increasing the complexity of the system from point unbalances to
distributed unbalance throughout the rotor as well as unbalances on the disks.
Finally in order to save power, the strategy of continuous use of AMBs has been changed to the use of
AMBs to get the required correction masses which when applied should prevent the unnecessary use of
power for continuous AMB use.
The FEM models have been solved using ODE23 in MATLAB and the results have been compiled for
various time ranges to show the vibration responses and attenuations with AMBs.
