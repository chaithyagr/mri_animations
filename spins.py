from manim import *
import numpy as np

N = 3
grid = np.asarray(np.meshgrid(np.arange(0, N), np.arange(0, N), np.arange(0, N),))

def get_spin(location, angle_in_xy, angle_in_z=PI/3, radius=0.2, color=RED, color_sphere=BLUE):
    sphere = Sphere(
                center=location,
                radius=radius,
                resolution=(16, 16),
                u_range=[0.001, PI - 0.001],
                v_range=[0, TAU],
                sheen_factor=0,
            )
    sphere.set_color(color_sphere)
    arrow_radius = radius * 4
    radius_in_plane = arrow_radius*np.cos(angle_in_z)
    radius_out_plane = arrow_radius*np.sin(angle_in_z)
    delta = [radius_in_plane*np.sin(angle_in_xy), radius_in_plane*np.cos(angle_in_xy), radius_out_plane]
    start_location = location - delta
    end_location = location + delta
    arrow = Arrow3D(
                start=start_location,
                end=end_location,
                resolution=8,
                color=color,
        )
    M0 = Group(sphere, arrow)
    return M0

class ThreeDSpinsRandomOrientation(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    arrow_angle = np.random.random()*2*PI
                    angle_z = np.random.random()*PI/2 
                    M0 = get_spin(grid[:, i, j, k], arrow_angle, angle_in_z=angle_z, radius=0.1)
                    self.add(axes, M0)
        self.renderer.camera.light_source.move_to(3*IN) # changes the source of the light
        self.set_camera_orientation(phi=75 * DEGREES, theta=30 * DEGREES)
        self.begin_3dillusion_camera_rotation(rate=2)
        self.wait()
        self.stop_3dillusion_camera_rotation()
        
        
class ThreeDSpinsPrecession(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        def rotate(d, dt):
            d.rotate(dt*5, about_point=d.get_center())
            
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    arrow_angle = np.random.random()*2*PI
                    M0 = get_spin(grid[:, i, j, k], arrow_angle, radius=0.1)
                    M0.add_updater(rotate)
                    self.add(axes, M0)
                    
        self.renderer.camera.light_source.move_to(3*IN) # changes the source of the light
        self.set_camera_orientation(phi=75 * DEGREES, theta=30 * DEGREES)
        self.wait(2*PI)        
        
        
class ParaSurface(ThreeDScene):
    def func(self, u, v):
        return np.array([np.cos(v)*np.sin(v), np.sin(u)*np.sin(v), u])

    def construct(self):
        axes = ThreeDAxes()
        surface = Surface(
            lambda u, v: axes.c2p(*self.func(u, v)),
            u_range=[-PI, PI],
            v_range=[0, TAU],
            resolution=8,
        )
        #surface.set_color(BLUE)
        self.set_camera_orientation(theta=70 * DEGREES, phi=75 * DEGREES)
        self.add(axes, surface)
        self.wait()
        self.set_camera_orientation(theta=15 * DEGREES, phi=15 * DEGREES)
