#Trial of standing waves

from standing import Standing
import pdb
import matplotlib.pyplot as plt

wave = Standing()

print(wave.constants)
print(wave.deltas)


sol = wave.solve()
plt.plot([n for n in range(len(sol))],[sol[i][0][0] for i in range(len(sol))])
plt.show()

print(wave.constants)
print(wave.deltas)
print(wave._js,wave.js)



