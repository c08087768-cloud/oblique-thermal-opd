"""Small solution interfaces composed with the existing numerical algorithms.

No verification, plotting, file writes or FEM work is triggered here.
All temperatures are rises above ambient; OPD follows the fixed refracted ray.
"""
from dataclasses import dataclass
from thermal_opd.steady import steady_temperature,steady_temperature_gradient,axial_response
from thermal_opd.opd import steady_refracted_opd,prepare_modal_ray_series
from thermal_opd.transient import TransientModel


def _ray_settings(prepared,num_ray_points,point_chunk_size):
    p = prepared.parameters
    return dict(radius=p.radius,refraction_angle=prepared.refraction_angle,
                num_ray_points=prepared.resolution.ray_points if num_ray_points is None else num_ray_points,
                point_chunk_size=point_chunk_size)


@dataclass(frozen=True)
class SteadySolution:
    prepared: object

    def _thermal_args(self):
        p = self.prepared.parameters
        return (self.prepared.coefficients,self.prepared.roots,p.conductivity,p.thickness,p.h_eff,p.h_eff)

    def temperature(self,r,phi,z):
        return steady_temperature(r,phi,z,*self._thermal_args())

    def temperature_gradient(self,r,phi,z):
        return steady_temperature_gradient(r,phi,z,*self._thermal_args())

    def opd(self,x,y,*,num_ray_points=None,point_chunk_size=512):
        return steady_refracted_opd(x,y,*self._thermal_args(),self.prepared.parameters.dn_dT,
                                    **_ray_settings(self.prepared,num_ray_points,point_chunk_size))

    def prepare_opd(self,*,num_ray_points=None,point_chunk_size=512):
        """Prepare shared axial weights once for repeated aperture queries."""
        p,c,roots = self.prepared.parameters,self.prepared.coefficients,self.prepared.roots
        evaluator = prepare_modal_ray_series(c,roots,
            lambda m,n,z:axial_response(z,roots[m][n],c[m][n],p.conductivity,p.thickness,p.h_eff,p.h_eff)[0][None,:],
            1,p.thickness,p.dn_dT,**_ray_settings(self.prepared,num_ray_points,point_chunk_size))
        return lambda x,y:evaluator(x,y)[0]


@dataclass(frozen=True,init=False)
class TransientSolution:
    prepared: object
    _solver: TransientModel

    def __init__(self,prepared,*,num_axial_modes,history=None):
        object.__setattr__(self,'prepared',prepared)
        object.__setattr__(self,'_solver',TransientModel(prepared.parameters,prepared.coefficients,prepared.roots,
                                                       num_axial_modes=num_axial_modes,history=history))

    @property
    def parameters(self):
        return self.prepared.parameters

    @property
    def coefficients(self):
        return self.prepared.coefficients

    @property
    def eigenvalues(self):
        return self.prepared.roots

    @property
    def history(self):
        return self._solver.history

    @property
    def axial(self):
        return self._solver.axial

    @property
    def slowest_decay_rate(self):
        return self._solver.slowest_decay_rate

    def temperature(self,r,phi,z,time):
        return self._solver.temperature(r,phi,z,time)

    def temperature_gradient(self,r,phi,z,time):
        return self._solver.temperature_gradient(r,phi,z,time)

    def opd(self,x,y,time,*,num_ray_points=None,point_chunk_size=512):
        return self._solver.refracted_opd(x,y,time,**_ray_settings(self.prepared,num_ray_points,point_chunk_size))

    def opd_series(self,x,y,times,*,num_ray_points=None,point_chunk_size=512):
        return self._solver.refracted_opd_series(x,y,times,**_ray_settings(self.prepared,num_ray_points,point_chunk_size))

    def prepare_opd(self,time,*,num_ray_points=None,point_chunk_size=512):
        return self._solver.prepare_refracted_opd(time,**_ray_settings(self.prepared,num_ray_points,point_chunk_size))
