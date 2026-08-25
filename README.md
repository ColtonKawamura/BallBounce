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


# Results with Appropriate Physical Parameters

With parameters that match those in the lab:

- $\hat{m}= 20$; the mass of the ball is much more than the mass of the paper-cells.
- $\hat{gamma} = [.0006, 0023]$; damping should be included, but is not a dominating factor
- $\hat{v}= .1$; impact kinetic energy is 10% of the potential energy of the springs.
- $\hat{g} = 0$; gravity is not a significant factor of the experiment

We get the following curves:

<img width="696" height="931" alt="image" src="https://github.com/user-attachments/assets/3f6a4664-25a7-4aa4-8e47-0e833925171f" />

The actual visualization of the simulations are shown in the following video.

- **White circles**  
  Default state for all particles: no special event has occurred yet.

- **Yellow circle (particle 2, just below the ball)**  
  Turns yellow when the **reflected wave is predicted to arrive back at the top** of the stack (based on the measured travel time to the bottom). This highlights the moment when energy carried by the wave reaches the upper contact again.

- **Red circle (bottom particle)**  
  Turns red once the **wave first reaches and compresses the bottom contact**. This marks the moment the initial impact pulse has traversed the entire stack.

- **Cyan circle (ball)**  
  The ball turns cyan after it has **lost contact with the stack** (gap larger than one diameter plus a small tolerance). This visually marks the onset of rebound / free flight.

https://github.com/user-attachments/assets/db07f7d6-c86a-4d45-a3be-91061ca162a0

# Pressure Effect

Decreasing pressure effectively decreases spring constant of the paper-air cell. Because effective damping and spring constant are coupled

$$
\hat{\gamma} = \frac{\gamma}{\sqrt{k_c\, m_c}}
$$

decreasing $k_c$ would also cause a decrease in $\hat{g}$. Plotting this explicit, we see,

<img width="693" height="926" alt="image" src="https://github.com/user-attachments/assets/cac70222-5204-4038-a43f-682acd66ee27" />

this causes a shift in peak to lower $N$ as well as decreasing overall restitution. 




