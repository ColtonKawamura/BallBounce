For the ball and paper sheet experiment, the anomalous height may be the result of a reflected wave from the initial impact.

This simulation provides an idealized model to help isolate this as the main effect.
Instead of paper sheets and a ball, we idealize the system as a 1D chain of particles
with damped harmonic contacts, impacted by a single free particle under gravity.

The following parameters can be tuned to investigate this effect:

- $\hat{m} = \frac{m_b}{m_c}$ =  ratio of ball to chain mass
- $\hat{k} = \frac{k_b}{k_c}$ = spring constant of the ball
- $\hat{d} = \frac{d_b}{d_c}$ diameter ratio of ball to chain
- $N$ = number of particles (already dimensionless)
- $\hat{\gamma} = \frac{\gamma}{\sqrt{k_c\, m_c}}$ 
	- When $\hat{\gamma} >1$ means overdamped springs.
- $\hat{v} = \frac{v_0}{d\sqrt{\frac{k_c}{m_c}}}$ =
	- When $\hat{v}>1$, means that the kinetic energy from impact is greater than potential energy of the spring.

$$
\hat{v} = \frac{v_0}{d\sqrt{\frac{k_c}{m_c}}} \rightarrow \hat{v}^2 = \frac{v_0^2\,m}{d^2k}
$$

- $\hat{g} = \frac{g}{d_c^2\frac{k_c}{m_c}}$ 
	- When $\hat{g}>1$, means that the gravity is stronger than force to full compress the chain spring

$$
\hat{g} = \frac{g}{d^2\frac{k_c}{m_c}} = \frac{gm}{d_ck_c}
$$

- $e$ = coefficient of restitution (dimensionless)
- $\tau$  = contact time (time -TBD)
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
