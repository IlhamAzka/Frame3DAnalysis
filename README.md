# Frame3DAnalysis

## Dividing Equation for LDTYPE 3
Because of the long equation of fixed end shear and fixed end moment on the beginning end of the member the equation is divided into different parts so the matlab code for the equation more readable

![equation](http://www.sciweavers.org/tex2img.php?eq=FS_b%20%3D%20%5Cunderbrace%7B%5Cfrac%7Bw_1%28L%20-%20l_1%29%5E3%7D%7B20L%5E3%7D%7D_%5Ctext%7BL%7D%5C%7B%5Cunderbrace%7B%287L%20%2B%208l_1%29%7D_%5Ctext%7BA%7D%5C%7D%20-%20%5Cunderbrace%7B%5Cfrac%7Bl_2%283L%20%2B%202l_1%29%7D%7BL%20-%20l_1%7D%7D_%5Ctext%7BB%7D%5B%5Cunderbrace%7B1%20%2B%20%5Cfrac%7Bl_2%7D%7BL%20-%20l_1%7D%20%2B%20%5Cfrac%7Bl_2%5E2%7D%7B%28L%20-%20l_1%29%5E2%7D%7D_%5Ctext%7BB1%7D%5D%20%2B%20%5Cunderbrace%7B%5Cfrac%7B2l_2%5E4%7D%7B%28L%20-%20l_1%29%5E3%7D%7D_%5Ctext%7BC%7D%5C%7D%20%2B%20%5Cunderbrace%7B%5Cfrac%7Bw_2%28L%20-%20l_1%29%5E3%7D%7B20L%5E3%7D%7D_%5Ctext%7BR%7D%5C%7B%5Cunderbrace%7B%283L%20%2B%202l_1%29%7D_%5Ctext%7BD%7D%5B%5Cunderbrace%7B1%20%2B%20%5Cfrac%7Bl_2%7D%7BL%20-%20l_1%7D%20%2B%20%5Cfrac%7Bl_2%5E2%7D%7B%28L%20-%20l_1%29%5E2%7D%7D_%5Ctext%7BD1%7D%5D%20-%20%5Cunderbrace%7B%5Cfrac%7Bl_2%5E3%7D%7B%28L%20-%20l_1%29%5E2%7D%7D_%5Ctext%7BE%7D%5Cunderbrace%7B2%20%2B%20%5Cfrac%7B15L%20-%208l_2%7D%7BL%20-%20l_1%7D%7D_%5Ctext%7BE1%7D%5C%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

![equation](http://www.sciweavers.org/tex2img.php?eq=FM_b%20%3D%20%5Cunderbrace%7B%5Cfrac%7Bw_1%28L%20-%20l_1%29%5E3%7D%7B60L%5E3%7D%7D_%5Ctext%7BL%7D%20%0A%5C%7B%20%5Cunderbrace%7B3%28L%20%2B%204l_1%29%7D_%5Ctext%7BA%7D%20-%20%0A%20%20%20%5Cunderbrace%7B%5Cfrac%7Bl_2%282L%20%2B%203l_1%29%7D%7BL%20-%20l_1%7D%7D_%5Ctext%7BB%7D%20%20%0A%5B%20%20%5Cunderbrace%7B1%20%2B%20%5Cfrac%7Bl_2%7D%7BL%20-%20l_1%7D%20%2B%20%5Cfrac%7Bl_2%5E2%7D%7B%28L%20-%20l_1%29%5E2%7D%7D_%5Ctext%7BB1%7D%20%5D%20%2B%0A%20%20%20%5Cunderbrace%7B%5Cfrac%7B3l_2%5E4%7D%7B%28L%20-%20l_1%29%5E3%7D%7D_%5Ctext%7BC%7D%20%0A%5C%7D%20%2B%0A%20%20%20%5Cunderbrace%7B%5Cfrac%7Bw_2%28L%20-%20l_1%29%5E3%7D%7B60L%5E3%7D%7D_%5Ctext%7BR%7D%0A%5C%7B%20%5Cunderbrace%7B%282L%20%2B%203l_1%29%7D_%5Ctext%7BD%7D%0A%5B%20%20%5Cunderbrace%7B1%20%2B%20%5Cfrac%7Bl_2%7D%7BL%20-%20l_1%7D%20%2B%20%5Cfrac%7Bl_2%5E2%7D%7B%28L%20-%20l_1%29%5E2%7D%7D_%5Ctext%7BD1%7D%20%5D%20-%0A%20%20%20%5Cunderbrace%7B%5Cfrac%7B3l_2%5E3%7D%7B%28L%20-%20l_1%29%5E2%7D%7D_%5Ctext%7BE%7D%0A%20%20%20%5Cunderbrace%7B1%20%2B%20%5Cfrac%7B5L%20-%204l_2%7D%7BL%20-%20l_1%7D%7D_%5Ctext%7BE1%7D%0A%5C%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

<!-- $$ -->
<!-- FS_b = \underbrace{\frac{w_1(L - l_1)^3}{20L^3}}_\text{L}  -->
<!-- \{ \underbrace{(7L + 8l_1)}_\text{A} \} -->
<!-- -  \underbrace{\frac{l_2(3L + 2l_1)}{L - l_1}}_\text{B}   -->
<!-- [  \underbrace{1 + \frac{l_2}{L - l_1} + \frac{l_2^2}{(L - l_1)^2}}_\text{B1} ] + -->
<!--    \underbrace{\frac{2l_2^4}{(L - l_1)^3}}_\text{C}  -->
<!-- \} + -->
<!--    \underbrace{\frac{w_2(L - l_1)^3}{20L^3}}_\text{R} -->
<!-- \{ \underbrace{(3L + 2l_1)}_\text{D} -->
<!-- [  \underbrace{1 + \frac{l_2}{L - l_1} + \frac{l_2^2}{(L - l_1)^2}}_\text{D1} ] - -->
<!--    \underbrace{\frac{l_2^3}{(L - l_1)^2}}_\text{E} -->
<!--    \underbrace{2 + \frac{15L - 8l_2}{L - l_1}}_\text{E1} -->
<!-- \} -->
<!-- $$ -->

<!-- $$FM_b =  -->
<!--    \underbrace{\frac{w_1(L - l_1)^3}{60L^3}}_\text{L}  -->
<!-- \{ \underbrace{3(L + 4l_1)}_\text{A} -  -->
<!--    \underbrace{\frac{l_2(2L + 3l_1)}{L - l_1}}_\text{B}   -->
<!-- [  \underbrace{1 + \frac{l_2}{L - l_1} + \frac{l_2^2}{(L - l_1)^2}}_\text{B1} ] + -->
<!--    \underbrace{\frac{3l_2^4}{(L - l_1)^3}}_\text{C}  -->
<!-- \} + -->
<!--    \underbrace{\frac{w_2(L - l_1)^3}{60L^3}}_\text{R} -->
<!-- \{ \underbrace{(2L + 3l_1)}_\text{D} -->
<!-- [  \underbrace{1 + \frac{l_2}{L - l_1} + \frac{l_2^2}{(L - l_1)^2}}_\text{D1} ] - -->
<!--    \underbrace{\frac{3l_2^3}{(L - l_1)^2}}_\text{E} -->
<!--    \underbrace{1 + \frac{5L - 4l_2}{L - l_1}}_\text{E1} -->
<!-- \}$$ -->
