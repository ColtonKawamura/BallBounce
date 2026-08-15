For the ball and paper sheet experiment, the anomalous height may be the result of a reflected wave from the initial impact.

This simulation provides an idealized model to help isolate this as the main effect.
Instead of paper sheets and a ball, we idealize the system as a 1D chain of particles
with damped harmonic contacts, impacted by a single free particle under gravity.

The following parameters can be tuned to investigate this effect:

- **`scalHeightDrop`** — drop height of the impacting particle, setting the impact velocity $v = \sqrt{2g \cdot h}$
- **`scalGravity`** — gravitational acceleration; controls how quickly particle 1 falls and returns after bouncing
- **`scalSpringConst`** — contact stiffness; stiffer contacts produce shorter, sharper wave pulses that reflect faster
- **`scalDamp`** — viscous damping; higher values dissipate the reflected wave before it can return to particle 1
- **`scalNumPart`** — number of particles in the chain; sets the distance the wave travels before reflecting off the fixed wall
- **`scalDiam`** — particle diameter; sets the lattice spacing and therefore the wave travel time
- **`scalMass`** — particle mass; together with `scalSpringConst` sets the natural frequency $\omega_0 = \sqrt{k/m}$ and wave speed
- **`scalPressure`** — static pre-compression; controls the equilibrium overlap and effective contact stiffness

# Idealized Effect

Below is an idealized effect. If the parameters are chosen in a certain way, the wave from the impact will travel through the chain and reflect back to eject the free particle upward.

<img width="1954" height="1174" alt="sim1d_construtve" src="https://github.com/user-attachments/assets/0940faa8-ed46-4947-a533-d3ed7f05b8f4" />

This results in a a significant rearwards (upward) velocity from the reflection paired with a recover of kinetic energy.

<img width="1527" height="919" alt="velocity" src="https://github.com/user-attachments/assets/013f1744-bb84-4c3d-a033-cc077b416786" />

# Non-Effect

If there are many layers to the stack, wave will either be dissipated by the time it would have reflect back or the wave will reach the impactor at a time that would result in constructive interference with it's oscillation from the impact.

Either way, the results in no significant recovery of kinetic energy.
<img width="1954" height="1174" alt="sim1d" src="https://github.com/user-attachments/assets/14d3381a-a7af-4117-bfdc-7a7fbdc6ee2a" />
<img width="1283" height="834" alt="nonconst_ke" src="https://github.com/user-attachments/assets/54813d03-68f1-42a1-b184-8b049be70885" />


This effect "knob" is expected to be a function of dissipation (damping factor, number of layers, etc). 

# Future Work

Intent is to run simulations and sweep over parameters to find a relation between increased kinetic energy recover and system parameters isolate main causes of the recovery.

- [ ] Add hertzian contacts
