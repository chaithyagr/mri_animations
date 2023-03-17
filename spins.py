from manim import *
#from manim.opengl import *
import numpy as np
from manim_slides import ThreeDSlide

N = 3
grid = np.asarray(np.meshgrid(np.arange(-N//2+1, N//2+1), np.arange(-N//2+1, N//2+1), np.arange(-N//2+1, N//2+1)))
RES = 2
np.random.seed(0)


def get_spin(location, angle_in_xy, angle_in_z=PI/3, radius=0.2, color=RED, sphere_color=BLUE):
    sphere = Sphere(
                center=location,
                radius=radius,
                resolution=(RES, RES),
                u_range=[0.001, PI - 0.001],
                v_range=[0, TAU],
            )
    sphere.set_color(sphere_color)
    arrow_radius = radius * 6
    radius_in_plane = arrow_radius*np.cos(angle_in_z)
    radius_out_plane = arrow_radius*np.sin(angle_in_z)
    delta = [radius_in_plane*np.sin(angle_in_xy), radius_in_plane*np.cos(angle_in_xy), radius_out_plane]
    start_location = location - delta
    end_location = location + delta
    arrow = Arrow3D(
                start=start_location,
                end=end_location,
                resolution=RES,
                color=color,
                thickness=0.001,
                height=0.1,
                base_radius=0.04,
        )
    return arrow, sphere

class ThreeDSpinsRandomOrientation(ThreeDSlide):
    def construct(self):
        axes = ThreeDAxes()
        arrows = np.empty((N, N, N), dtype=object)
        spheres =  np.empty((N, N, N), dtype=object)
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    arrow_angle = np.random.random()*2*PI
                    angle_z = np.random.random()*PI/2 
                    arrows[i, j, k], spheres[i, j, k] = get_spin(grid[:, i, j, k], arrow_angle, angle_in_z=angle_z, radius=0.5)
        self.add(*arrows, *spheres)
        self.renderer.camera.light_source.move_to(3*IN) # changes the source of the light
        self.set_camera_orientation(phi=75 * DEGREES, theta=30 * DEGREES)
        self.begin_3dillusion_camera_rotation(rate=2)
        self.wait(PI*2)
        self.stop_3dillusion_camera_rotation()
        self.next_slide()
        
class ThreeDSpinsPrecession(ThreeDSlide):
    def do_camera_rotation(self):
        self.begin_ambient_camera_rotation(rate=DEGREES*5, about='theta')
        self.wait(DEGREES*180)
        self.stop_ambient_camera_rotation()
        self.begin_ambient_camera_rotation(rate=-DEGREES*5, about='phi')
        self.wait(DEGREES*180)
        self.stop_ambient_camera_rotation()
        
        self.begin_ambient_camera_rotation(rate=-DEGREES*5, about='theta')
        self.wait(DEGREES*180)
        self.stop_ambient_camera_rotation()
        self.begin_ambient_camera_rotation(rate=DEGREES*5, about='phi')
        self.wait(DEGREES*180)
        self.stop_ambient_camera_rotation()
        
        
    def construct(self):
        axes = ThreeDAxes()
        self.add(axes)
        def rotate(d, dt, about='center'):
            if about == 'center':
                d.rotate(dt*5, about_point=d.get_center())
            elif about == 'start':
                d.rotate(dt*5, about_point=d.get_start())
        
        arrows = np.empty((N, N, N), dtype=object)
        spheres =  np.empty((N, N, N), dtype=object)
        Animations = []
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    arrow_angle = np.random.random()*2*PI
                    arrows[i, j, k], spheres[i, j, k] = get_spin(grid[:, i, j, k], arrow_angle, radius=0.1)
                    arrows[i, j, k].add_updater(rotate)
        self.add(*arrows.flatten(), *spheres.flatten())
        self.renderer.camera.light_source.move_to(3*IN) # changes the source of the light
        self.set_camera_orientation(phi=75 * DEGREES, theta=30 * DEGREES)
        self.start_loop()
        #self.do_camera_rotation() 
        #self.end_loop()  
        
        self.remove(*spheres.flatten())
        Animations = []
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    arrow = arrows[i, j, k]
                    line = Line(start=arrow.get_start(), end=[0,0,0])
                    Animations.append(MoveAlongPath(arrow, line))
        self.play(*Animations, rate_func=linear, run_time=2)
        
        self.move_camera(phi=80*DEGREES, theta=45 * DEGREES, zoom=3, run_time=2)
        self.start_loop() 
        self.wait(PI)
        self.end_loop()                    
        
        



