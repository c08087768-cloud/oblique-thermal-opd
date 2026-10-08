"""Independent 3D P1-triangle x P1-line thermal FEM on a polygonal cylinder.

Six-node prismatic trial functions: each 3D node has an independent unknown.
Constant k and face-specific Robin h permit exact tensor-product assembly.
The Gaussian is integrated directly on the front, without modal projection.
"""

from dataclasses import dataclass, field
from operator import index
from time import perf_counter

import numpy as np
from matplotlib.tri import Triangulation
from scipy.sparse import diags, kron
from scipy.sparse.linalg import cg, spsolve, factorized
from skfem import MeshTri, MeshLine, Basis, FacetBasis, ElementTriP1, ElementLineP1, LinearForm, BilinearForm, asm
from skfem.models.poisson import laplace


@BilinearForm
def value_mass(u, v, w):
    # scikit-fem 8's DiscreteField subclass predates NumPy 2's array_wrap API.
    return np.asarray(u)*np.asarray(v)


@dataclass
class FemSolution:
    disk_mesh: MeshTri
    z_m: np.ndarray
    temperature_K: np.ndarray  # (z node, xy node), rise above ambient
    diagnostics: dict
    _finder: object = field(init=False, repr=False)
    _origins: np.ndarray = field(init=False, repr=False)
    _inverse: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        p, t = self.disk_mesh.p, self.disk_mesh.t
        self._finder = Triangulation(p[0], p[1], t.T).get_trifinder()
        triangles = p[:, t].transpose(2, 1, 0)
        self._origins = triangles[:, 0]
        self._inverse = np.linalg.inv(np.stack((triangles[:, 1]-triangles[:, 0],
                                                triangles[:, 2]-triangles[:, 0]), axis=2))

    def transmitted_mask(self, x0, y0, refraction_angle):
        if not np.isfinite(refraction_angle) or abs(refraction_angle) >= np.pi/2:
            raise ValueError('Require a finite forward refraction angle')
        x, y = np.broadcast_arrays(x0, y0)
        shift = self.z_m[-1]*np.tan(refraction_angle)
        # Convex polygon: both endpoints inside guarantees the whole ray inside.
        return (self._finder(x, y) >= 0) & (self._finder(x, y+shift) >= 0)

    def temperature_at(self, x, y, z):
        """Interpolate the actual prism shape functions; reject extrapolation."""
        x, y, z = np.broadcast_arrays(x, y, z)
        if not all(np.all(np.isfinite(a)) for a in (x, y, z)):
            raise ValueError('Coordinates must be finite')
        cells = np.asarray(self._finder(x.ravel(), y.ravel()), dtype=int)
        if np.any(cells < 0) or np.any(z < 0) or np.any(z > self.z_m[-1]):
            raise ValueError('Point is outside the FEM polygon or axial interval')
        planar = np.column_stack((x.ravel(), y.ravel()))-self._origins[cells]
        bary = np.einsum('nij,nj->ni', self._inverse[cells], planar)
        bary = np.column_stack((1-bary.sum(axis=1), bary))
        nodes = self.disk_mesh.t[:, cells].T
        layer = np.clip(np.searchsorted(self.z_m, z.ravel(), side='right')-1, 0, self.z_m.size-2)
        fraction = (z.ravel()-self.z_m[layer])/(self.z_m[layer+1]-self.z_m[layer])
        front = np.sum(bary*self.temperature_K[layer[:, None], nodes], axis=1)
        rear = np.sum(bary*self.temperature_K[(layer+1)[:, None], nodes], axis=1)
        return ((1-fraction)*front+fraction*rear).reshape(x.shape)

    def refracted_opd(self, x0, y0, refraction_angle, dn_dT, *, point_chunk_size):
        """Integrate prism temperature exactly along each fixed straight ray.

        Split at all crossed xy edges and z layers. Within each segment the
        triangle-linear x line-linear field is quadratic in z; two Gauss
        points therefore integrate it exactly, up to floating point error.
        """
        chunk = index(point_chunk_size)
        if chunk < 1:
            raise ValueError('point_chunk_size must be positive')
        x, y = np.broadcast_arrays(x0, y0)
        if not np.all(self.transmitted_mask(x, y, refraction_angle)):
            raise ValueError('A ray leaves the FEM mesh before the rear face')
        slope = np.tan(refraction_angle)
        edges = self.disk_mesh.p[:, self.disk_mesh.facets]
        dx = edges[0, 1]-edges[0, 0]
        nonvertical = np.abs(dx) > 8*np.finfo(float).eps*np.max(np.abs(self.disk_mesh.p))
        xa, ya, delta_x, delta_y = (edges[0, 0, nonvertical], edges[1, 0, nonvertical],
                                    dx[nonvertical], (edges[1, 1]-edges[1, 0])[nonvertical])
        result = np.empty(x.size)
        for begin in range(0, x.size, chunk):
            all_z, all_x, all_y, all_weights, all_ids = [], [], [], [], []
            for j in range(begin, min(begin+chunk, x.size)):
                xx, yy = x.ravel()[j], y.ravel()[j]
                cuts = self.z_m
                if slope != 0:
                    fraction = (xx-xa)/delta_x
                    active = (fraction >= 0) & (fraction <= 1)
                    crossings = (ya[active]+fraction[active]*delta_y[active]-yy)/slope
                    crossings = crossings[(crossings > 0) & (crossings < self.z_m[-1])]
                    cuts = np.unique(np.concatenate((cuts, crossings)))
                middle = (cuts[1:]+cuts[:-1])/2
                half = (cuts[1:]-cuts[:-1])/2
                zq = (middle[:, None] + half[:, None]*np.array([-1., 1.])/np.sqrt(3)).ravel()
                all_z.append(zq)
                all_x.append(np.full(zq.size, xx))
                all_y.append(yy+slope*zq)
                all_weights.append(np.repeat(half, 2))
                all_ids.append(np.full(zq.size, j-begin, dtype=int))
            if not all_z:
                continue
            values = self.temperature_at(np.concatenate(all_x), np.concatenate(all_y), np.concatenate(all_z))
            totals = np.bincount(np.concatenate(all_ids), weights=values*np.concatenate(all_weights),
                                 minlength=min(chunk, x.size-begin))
            result[begin:begin+totals.size] = dn_dT*totals/np.cos(refraction_angle)
        return result.reshape(x.shape)


@dataclass
class FemSystem:
    disk_mesh: MeshTri
    z_m: np.ndarray
    matrix: object
    mass: object  # geometric consistent mass; multiply by rho*cp for heat capacity
    rhs: np.ndarray
    Mxy: object
    Mz: object
    Sxy: object
    diagnostics: dict


def assemble_thermal_fem(radius, thickness, conductivity, h_front, h_back, h_side, front_flux,
                         *, disk_refinements, axial_layers, quadrature_order):
    """Shared steady/transient prism space, conduction, Robin and front load."""
    if not all(np.isfinite(v) and v > 0 for v in (radius, thickness, conductivity)):
        raise ValueError('Geometry and conductivity must be finite and positive')
    if not all(np.isfinite(v) and v >= 0 for v in (h_front, h_back, h_side)) or max(h_front, h_back, h_side)==0:
        raise ValueError('Robin coefficients must be nonnegative, with some heat loss')
    level, layers, order = index(disk_refinements), index(axial_layers), index(quadrature_order)
    if min(level, layers) < 1 or order < 2:
        raise ValueError('Require refinements/layers >= 1 and quadrature >= 2')
    started = perf_counter()
    mesh = MeshTri.init_circle(nrefs=level).scaled(radius)
    z = np.linspace(0, thickness, layers+1)
    xy_basis = Basis(mesh, ElementTriP1(), intorder=order)
    side_basis = FacetBasis(mesh, ElementTriP1(), intorder=order)
    z_basis = Basis(MeshLine(z[None, :]), ElementLineP1(), intorder=order)
    Mxy, Kxy, Sxy = asm(value_mass, xy_basis), asm(laplace, xy_basis), asm(value_mass, side_basis)
    Mz, Kz = asm(value_mass, z_basis), asm(laplace, z_basis)
    first, last = np.zeros(layers+1), np.zeros(layers+1)
    first[0], last[-1] = 1., 1.
    matrix = (conductivity*kron(Kz, Mxy) + kron(Mz, conductivity*Kxy+h_side*Sxy)
              + kron(diags(h_front*first+h_back*last), Mxy)).tocsr()
    @LinearForm
    def source(v, w):
        q = front_flux(w.x[0], w.x[1])
        if not np.all(np.isfinite(q)):
            raise ValueError('Front heat flux must be finite')
        return q*np.asarray(v)
    planar_rhs = asm(source, xy_basis)
    rhs = np.zeros(matrix.shape[0])
    rhs[:mesh.nvertices] = planar_rhs
    assembly_seconds = perf_counter()-started
    edges = mesh.p[:, mesh.facets[:, mesh.boundary_facets()]]
    area, perimeter = float(np.sum(xy_basis.dx)), float(np.sum(side_basis.dx))
    inradius = float(np.min(np.linalg.norm((edges[:, 0]+edges[:, 1])/2, axis=0)))
    metadata = {'degrees_of_freedom': int(matrix.shape[0]), 'disk_triangles': int(mesh.nelements),
                'prism_elements': int(mesh.nelements*layers), 'disk_refinements': level,
                'axial_layers': layers, 'quadrature_order': order, 'assembly_seconds': assembly_seconds,
                'polygon_area_m2': area, 'polygon_perimeter_m': perimeter, 'polygon_inradius_m': inradius,
                'relative_area_deficit': float(1-area/(np.pi*radius**2)),
                'relative_perimeter_deficit': float(1-perimeter/(2*np.pi*radius))}
    return FemSystem(mesh,z,matrix,kron(Mz,Mxy).tocsr(),rhs,Mxy,Mz,Sxy,metadata)


def solve_steady_fem(radius, thickness, conductivity, h_front, h_back, h_side, front_flux,
                     *, disk_refinements, axial_layers, quadrature_order, solver, rtol):
    """Solve k*grad + Robin = prescribed front surface flux (SI)."""
    if solver not in ('direct','cg') or not 0 < rtol < 1:
        raise ValueError('Require direct/cg solver and 0 < rtol < 1')
    system = assemble_thermal_fem(radius,thickness,conductivity,h_front,h_back,h_side,front_flux,
                                  disk_refinements=disk_refinements,axial_layers=axial_layers,
                                  quadrature_order=quadrature_order)
    mesh, z, matrix, rhs = system.disk_mesh, system.z_m, system.matrix, system.rhs
    Mxy, Mz, Sxy = system.Mxy, system.Mz, system.Sxy
    solving = perf_counter()
    iterations = 0
    if solver == 'direct':
        vector = spsolve(matrix, rhs)
    else:
        def callback(_):
            nonlocal iterations
            iterations += 1
        vector, info = cg(matrix, rhs, M=diags(1/matrix.diagonal()), rtol=rtol, atol=0,
                          maxiter=10000, callback=callback)
        if info != 0:
            raise RuntimeError(f'FEM CG failed to converge: info={info}')
    solve_seconds = perf_counter()-solving
    if not np.all(np.isfinite(vector)):
        raise RuntimeError('FEM solution contains nonfinite values')
    temperature = vector.reshape(z.size, mesh.nvertices)
    front_loss = h_front*np.sum(Mxy@temperature[0])
    rear_loss = h_back*np.sum(Mxy@temperature[-1])
    side_loss = h_side*np.asarray(Mz.sum(axis=0)).ravel() @ temperature @ np.asarray(Sxy.sum(axis=0)).ravel()
    power = float(rhs.sum())
    residual = np.linalg.norm(matrix@vector-rhs)
    # Signed sources can have zero net input; preserve absolute diagnostics too.
    scale = max(abs(power), np.finfo(float).tiny)
    diagnostics = {**system.diagnostics, 'solver': solver, 'solver_rtol': rtol, 'iterations': iterations,
                   'solve_seconds': solve_seconds,
                   'input_W': power, 'front_loss_W': float(front_loss), 'rear_loss_W': float(rear_loss),
                   'side_loss_W': float(side_loss), 'balance_absolute_error_W': float(abs(front_loss+rear_loss+side_loss-power)),
                   'relative_energy_imbalance': float(abs(front_loss+rear_loss+side_loss-power)/scale),
                   'linear_residual_norm': float(residual),
                   'relative_linear_residual': float(residual/max(np.linalg.norm(rhs), np.finfo(float).tiny))}
    return FemSolution(mesh, z, temperature, diagnostics)


@dataclass
class TransientFemSolution:
    times_s: np.ndarray
    frames: list
    steps: list
    diagnostics: dict


def solve_transient_fem(radius, thickness, conductivity, h_front, h_back, h_side, front_flux,
                        *, density, specific_heat, times, time_step, history,
                        disk_refinements, axial_layers, quadrature_order, solver, rtol):
    """Zero initial rise, consistent capacity, backward Euler at exact observations.

    Each interval ends at the next requested observation (including fractional
    final steps).  Endpoint g(t_next) drives the implicit step.  Sparse factors
    are reused by dt.  Frames reuse FemSolution's interpolation/ray integration.
    """
    times = np.asarray(times,dtype=float)
    if times.ndim!=1 or not times.size or not np.all(np.isfinite(times)) or times[0]<0 or np.any(np.diff(times)<=0):
        raise ValueError('Observation times must be finite, nonnegative and strictly increasing')
    if not all(np.isfinite(v) and v>0 for v in (density,specific_heat,time_step)):
        raise ValueError('Capacity and time step must be finite and positive')
    if solver not in ('direct','cg') or not 0 < rtol < 1:
        raise ValueError('Require direct/cg solver and 0 < rtol < 1')
    system = assemble_thermal_fem(radius,thickness,conductivity,h_front,h_back,h_side,front_flux,
                                  disk_refinements=disk_refinements,axial_layers=axial_layers,
                                  quadrature_order=quadrature_order)
    capacity = density*specific_heat*system.mass
    vector = np.zeros_like(system.rhs)
    factors, steps, frames = {}, [], []
    time, power = 0., float(system.rhs.sum())
    started = perf_counter()
    for target in times:
        # Use an equal partition within each observation interval; avoids tiny
        # floating-point remainder steps and never exceeds the requested dt.
        count = int(np.ceil(max(0.,target-time)/time_step))
        if count:
            dt = (target-time)/count
            start_time = time
            step_matrix = (capacity/dt+system.matrix).tocsr()
            if solver=='direct' and dt not in factors:
                factors[dt] = factorized(step_matrix.tocsc())
            preconditioner = diags(1/step_matrix.diagonal()) if solver=='cg' else None
            for j in range(count):
                time = float(start_time+(j+1)*dt) if j+1<count else float(target)
                load = system.rhs*history.value(time)
                rhs = capacity@vector/dt+load
                previous = vector
                if solver=='direct':
                    vector = factors[dt](rhs)
                else:
                    vector, info = cg(step_matrix,rhs,x0=previous,M=preconditioner,rtol=rtol,atol=0.,maxiter=10000)
                    if info!=0:
                        raise RuntimeError(f'Transient FEM CG failed: info={info}')
                if not np.all(np.isfinite(vector)):
                    raise RuntimeError('Transient FEM has nonfinite temperatures')
                storage = float(np.sum(capacity@(vector-previous))/dt)
                loss = float(np.sum(system.matrix@vector))
                incoming = float(load.sum())
                residual = float(np.linalg.norm(step_matrix@vector-rhs)/max(np.linalg.norm(rhs),np.finfo(float).tiny))
                steps.append({'time_s':time,'dt_s':float(dt),'input_W':incoming,'storage_rate_W':storage,
                              'heat_loss_W':loss,'stored_energy_J':float(np.sum(capacity@vector)),
                              'balance_absolute_error_W':abs(storage+loss-incoming),
                              'relative_linear_residual':residual})
        diagnostic = {'time_s':float(target),'stored_energy_J':float(np.sum(capacity@vector))}
        frame = FemSolution(system.disk_mesh,system.z_m,vector.reshape(system.z_m.size,-1).copy(),diagnostic)
        if frames:
            # Geometry is identical; reuse the first frame's spatial locator.
            frame._finder,frame._origins,frame._inverse = frames[0]._finder,frames[0]._origins,frames[0]._inverse
        frames.append(frame)
    return TransientFemSolution(times,frames,steps,{**system.diagnostics,'time_integrator':'backward_euler',
        'time_step_max_s':float(time_step),'solver':solver,'solver_rtol':rtol,'density_kg_per_m3':density,
        'specific_heat_J_per_kg_K':specific_heat,'full_input_W':power,'solve_seconds':perf_counter()-started})
