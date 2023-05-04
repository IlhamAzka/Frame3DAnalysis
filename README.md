# Frame3DAnalysis

## Dividing Equation for LDTYPE 4
Because of the long equation of fixed end shear and fixed end moment on the beginning end of the member the equation is divided into different parts so the matlab code for the equation more readable

equation:

<img src="fsb_eq_black_bg.svg">

<img src="fmb_eq_black_bg.svg">

<!-- $$
FM_b = 
   \underbrace{\frac{w_1(L - l_1)^3}{20L^3}}_\text{L} 
\{ \underbrace{(7L + 8l_1)}_\text{A} - 
   \underbrace{\frac{l_2(3L + 2l_1)}{L - l_1}}_\text{B}  
[  \underbrace{1 + \frac{l_2}{L - l_1} + \frac{l_2^2}{(L - l_1)^2}}_\text{B1} ] +
   \underbrace{\frac{2l_2^4}{(L - l_1)^3}}_\text{C} 
\} +
   \underbrace{\frac{w_2(L - l_1)^3}{20L^3}}_\text{R}
\{ \underbrace{(3L + 2l_1)}_\text{D}
[  \underbrace{1 + \frac{l_2}{L - l_1} + \frac{l_2^2}{(L - l_1)^2}}_\text{D1} ] -
   \underbrace{\frac{l_2^3}{(L - l_1)^2}}_\text{E}
[  \underbrace{1 + \frac{15L - 8l_2}{L - l_1}}_\text{E1} ]
\}
$$

$$
FM_b = 
   \underbrace{\frac{w_1(L - l_1)^3}{60L^3}}_\text{L} 
\{ \underbrace{3(L + 4l_1)}_\text{A} - 
   \underbrace{\frac{l_2(2L + 3l_1)}{L - l_1}}_\text{B}  
[  \underbrace{1 + \frac{l_2}{L - l_1} + \frac{l_2^2}{(L - l_1)^2}}_\text{B1} ] +
   \underbrace{\frac{3l_2^4}{(L - l_1)^3}}_\text{C} 
\} +
   \underbrace{\frac{w_2(L - l_1)^3}{60L^3}}_\text{R}
\{ \underbrace{(2L + 3l_1)}_\text{D}
[  \underbrace{1 + \frac{l_2}{L - l_1} + \frac{l_2^2}{(L - l_1)^2}}_\text{D1} ] -
   \underbrace{\frac{3l_2^3}{(L - l_1)^2}}_\text{E}
[  \underbrace{1 + \frac{5L - 4l_2}{L - l_1}}_\text{E1} ]
\}
$$ -->
