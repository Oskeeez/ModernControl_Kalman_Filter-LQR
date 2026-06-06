# Mag-lev Train Control System

LQR control system for a two-carriage magnetic levitation train, modelled and simulated
in MATLAB. The controller reaches and maintains cruise speed while keeping passenger
acceleration within comfort limits and preventing the carriages from diverging from their
coupled positions.

## System

Two carriages coupled by a nonlinear spring-damper assembly. The carriages have
different masses (9,000 kg and 10,000 kg) and different aerodynamic drag profiles,
making symmetric control impossible - the controller has to account for the asymmetry
directly.

| Parameter | Carriage 1 | Carriage 2 |
|-----------|-----------|-----------|
| Mass | 9,000 kg | 10,000 kg |
| Drag coefficient | 0.5 N s²/m² | 0.45 N s²/m² |

Coupler: linear spring 5×10⁴ N/m, linear damping 1×10⁵ N s/m, nonlinear spring
1×10⁴ N/m³, natural length 10 m.

## Design

**Modelling**
The full system was modelled in state space with states for carriage positions and
velocities. The coupler force includes a cubic nonlinear term to capture the stiffening
behaviour at large displacements.

**Controllability problem**
After linearising about the steady-state operating point and computing the Jacobian,
the controllability matrix came back rank 3 of 4 - the system as formulated was not
fully controllable. The fix was to reduce the state representation by expressing carriage
separation as a single state rather than tracking both positions independently. The
reduced system achieved full rank and the design could proceed.

**LQR with integral action**
An LQR controller was designed on the reduced system and tuned to balance speed of
response against actuator effort. Integral action was added to eliminate steady-state
error - without it the train settles at a speed offset from the target due to the
asymmetric drag.

**Output clamping and anti-windup**
The controller output is clamped to enforce a maximum acceleration limit for passenger
comfort. When the output saturates, a standard anti-windup scheme prevents the
integrator from accumulating error it cannot act on, which would otherwise cause
significant overshoot when the saturation clears.

**Observer**
A state observer was designed to estimate the full state vector from the available
velocity outputs, allowing the LQR feedback law to operate without direct position
measurement.

## Requirements

- MATLAB R2020a or later
- Control System Toolbox
