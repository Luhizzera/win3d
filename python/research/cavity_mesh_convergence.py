"""Convergencia de malha da cavidade nao-estruturada (DIVIDA_TECNICA 5.3)."""
import sys, math, time
sys.path.insert(0,'build/bindings/python/Release')
import aether_core_py as core, aether_mesh_py as m, aether_solver_py as sol
def suite_lattice(n):
    t=m.DelaunayTetrahedralization3D(); step=1.0/n; jit=0.22*step
    for i in range(n+1):
        for j in range(n+1):
            for k in range(n+1):
                x,y,z=i*step,j*step,k*step
                if 0<i<n: x+=jit*math.sin(4.1*y+2.7*z+1.3)
                if 0<j<n: y+=jit*math.sin(3.3*z+5.1*x+0.7)
                if 0<k<n: z+=jit*math.sin(2.9*x+4.7*y+2.1)
                t.add_point(x,y,z)
    t.tetrahedralize(); return m.TetrahedralMesh.from_tetrahedralization(t)
PI=math.pi
def lid(p):
    if p.z < 1.0-1e-9: return core.Vector3(0,0,0)
    return core.Vector3(math.sin(PI*p.x)**2 * math.sin(PI*p.y)**2, 0.0, 0.0)
# Quantidades integrais, nao pontuais: media ponderada por volume nas faixas
# junto a tampa e ao fundo, e energia cinetica total. Sao o que a topologia
# afirma e o que converge com a malha.
def measure(mesh,s):
    tn=td=bn=bd=ke=0.0
    for c in range(mesh.cell_count()):
        v=mesh.cell_volume(c); u=s.velocity(c)
        ke += v*(u.x*u.x+u.y*u.y+u.z*u.z)
        z=mesh.cell_centroid(c).z
        if z>0.8: tn+=v*u.x; td+=v
        elif z<0.2: bn+=v*u.x; bd+=v
    return tn/td, bn/bd, 0.5*ke
print(' n  celulas   dt        passos   u topo     u fundo    E cinetica  div       t(s)')
prev=None
for n in (3,4,5,6,7):
    t0=time.time(); mesh=suite_lattice(n); tmesh=time.time()-t0
    s=sol.UnstructuredCavitySolver3D(mesh,0.1,lid)
    dt=s.stable_time_step(); steps=int(8.0/dt)
    t0=time.time()
    for _ in range(steps): s.step(dt)
    el=time.time()-t0
    top,bot,ke=measure(mesh,s)
    print('%2d %7d  %.2e %7d  %+.6f  %+.6f  %.6e  %.2e  %.0f (+%.0f malha)'
          % (n,mesh.cell_count(),dt,steps,top,bot,ke,s.max_face_divergence(),el,tmesh))
    if prev:
        print('     variacao vs anterior: topo %+.2f%%  fundo %+.2f%%  E %+.2f%%'
              % (100*(top-prev[0])/abs(prev[0]), 100*(bot-prev[1])/abs(prev[1]),
                 100*(ke-prev[2])/prev[2]))
    prev=(top,bot,ke)
