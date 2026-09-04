import matplotlib.pyplot as plt
import numpy as np 

t  = np.arange(20)
u = t / t[-1]  # or t / 19

# Parabolic decay from 20 to 0.1
#n = 20 - (20 - 0.1) * (u**2)
#p = 2 - (2 - 0.01) * (u**2)
#si = 10 - (10 - 0.05) * (u**2)

n  = 15 - 14.9 * (t / (t + 1)) 
p  = 1.5 - 1.49 * (t / (t + 0.025)) 
si = 10 - 9.95 * (t / (t + 0.025)) 

#n = n*0+0.1
#p = p*0+0.001
#si = si*0+0.05

plt.figure()
ax1 = plt.subplot(2,1,1)
ax1.plot(t,n,label='Nitrate')
ax1.plot(t,p,label='Phosphate')
ax1.plot(t,si,label='Silicate')
ax1.legend()
#plt.show()

nlimdia = np.zeros(20)
nlimfla = np.zeros(20)
nlimcoc = np.zeros(20)

rNdia = 0.5  ; rNfla = 0.5 ; rNcoc = 1.0
rPdia = 0.05 ; rPfla = 0.05 ; rPcoc = 0.0015
rSidia = 0.5 ; rSifla = 0.0 ; rSicoc = 0.0

upNdia = n/(n + rNdia)
upPdia = p/(p + rPdia)
upSidia = si/(si + rSidia)
upNfla = n/(n + rNfla)
upPfla = p/(p + rPfla)
upSifla = np.ones(20) # No silicate requirement for flagellates
upNcoc = n/(n + rNcoc)
upPcoc = p/(p + rPcoc)
upSicoc = np.ones(20) # No silicate requirement for coccolithophores

for i in range(20):
    nlimdia[i] = min(upNdia[i],upPdia[i],upSidia[i])
    nlimfla[i] = min(upNfla[i],upPfla[i],upSifla[i])
    nlimcoc[i] = min(upNcoc[i],upPcoc[i],upSicoc[i])

mu_dia = 0.62   ; mu_fla = 0.43   ; mu_coc = 0.4
#mu_dia = 0.8   ; mu_fla = 0.6   ; mu_coc = 0.74
TdepCoc = 0.4 * 0.1419 * t**0.8151
TdepFla = 0.43 * 1.67**(t/10.0)
TdepDia = 0.62 * 1.55**(t/10.0)    

grow_dia = TdepDia * mu_dia * nlimdia
grow_fla = TdepFla * mu_fla * nlimfla
grow_coc = TdepCoc * mu_coc * nlimcoc


ax2 = plt.subplot(2,1,2)
ax2.plot(t,grow_dia,label='Diatom')
ax2.plot(t,grow_fla,label='Flagellate')
ax2.plot(t,grow_coc,label='Coccolithophore')
ax2.legend()
plt.show()