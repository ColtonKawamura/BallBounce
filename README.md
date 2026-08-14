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

<img width="1527" height="919" alt="kinetic_energy" src="https://github.com/user-attachments/assets/b375e364-eb35-44de-a367-74799eaeede6" />
<img width="1527" height="919" alt="velocity" src="https://github.com/user-attachments/assets/013f1744-bb84-4c3d-a033-cc077b416786" />

# Damped Effect

If there are many layers to the stack, wave will be dissipated by the time it would have reflect back, resulting in no recovery of kinetic energy.
<img width="1954" height="1174" alt="sim1d_longchain" src="https://github.com/user-attachments/assets/249fcb75-14e7-4ad3-aae3-560978aa4add" />
<img width="1527" height="919" alt="velocity_longChain" src="https://github.com/user-attachments/assets/d348cb83-1804-4090-b005-6ad0236d6e4b" />
<img width="1527" height="919" alt="kinetic_energy_longChain" src="https://github.com/user-attachments/assets/55583794-8961-4b1c-9270-8ed0577bdcaf" />


This effect "knob" is expected to be a function of dissipation (damping factor, number of layers, etc). 

# Future Work

Intent is to run simulations and sweep over parameters to find a relation between increased kinetic energy recover and system parameters isolate main causes of the recovery.
