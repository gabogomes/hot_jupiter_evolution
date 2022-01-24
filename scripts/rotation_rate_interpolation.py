import numpy as np
import matplotlib.pyplot as plt

x=[0.5 , 0.8 , 1.0]
y=[3.0 , 8.0 , 14.0]
pol=np.polyfit(x,y,2)

xx = np.linspace(min(x),max(x))
yy = np.polyval(pol,xx)        

print(pol)

d_x=0.1
d_y=2.0

plt.plot(xx, yy, '-',x, y, 'ro')
plt.axis([min(xx)-d_x, max(xx)+d_x, min(yy)-d_y, max(yy)+d_y])
plt.xlabel(r"$M_{\star}$ [$M_{\odot}$]")
plt.ylabel(r"$\omega_{sat}$ [$\Omega_{\odot}$]")
plt.tight_layout()
plt.savefig("interpolation_rotation_vs_mass.pdf")
plt.show()