from manim import *
#from manim.opengl import *
import numpy as np
from manim_slides import ThreeDSlide

N = 3
grid = np.asarray(np.meshgrid(np.arange(-N//2+1, N//2+1), np.arange(-N//2+1, N//2+1), np.arange(-N//2+1, N//2+1)))
RES = 8
np.random.seed(0)
TEST = True

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


def get_main_spin():
    M0 = Arrow3D(
            start=[0, 0, 0],
            end=[0, 0, 1],
            resolution=RES,
            color=WHITE,
            thickness=0.05,
            height=0.2,
            base_radius=0.1,
        )
    return M0

def do_camera_rotation(obj):
    if TEST == False:
        obj.begin_ambient_camera_rotation(rate=DEGREES*5, about='theta')
        obj.wait(DEGREES*180)
        obj.stop_ambient_camera_rotation()
        obj.begin_ambient_camera_rotation(rate=-DEGREES*5, about='phi')
        obj.wait(DEGREES*180)
        obj.stop_ambient_camera_rotation()
        
        obj.begin_ambient_camera_rotation(rate=-DEGREES*5, about='theta')
        obj.wait(DEGREES*180)
        obj.stop_ambient_camera_rotation()
        obj.begin_ambient_camera_rotation(rate=DEGREES*5, about='phi')
        obj.wait(DEGREES*180)
        obj.stop_ambient_camera_rotation()
    else:
        obj.wait()


def add_axes(obj):
    axes = ThreeDAxes()
    obj.add(axes)
    obj.renderer.camera.light_source.move_to(3*IN) # changes the source of the light
    obj.set_camera_orientation(phi=75 * DEGREES, theta=30 * DEGREES)


def setup_all_spins(obj, updater=None, locations='grid'):
    arrows = np.empty((N, N, N), dtype=object)
    spheres =  np.empty((N, N, N), dtype=object)
    if locations == 'grid':
        loc = grid
    else:
        loc = locations
    for i in range(N):
        for j in range(N):
            for k in range(N):
                arrow_angle = np.random.random()*2*PI
                arrows[i, j, k], spheres[i, j, k] = get_spin(loc[:, i, j, k], arrow_angle, radius=0.1)
                if updater is not None:
                    arrows[i, j, k].add_updater(updater)
    obj.add(*arrows.flatten(), *spheres.flatten())
    return arrows, spheres

def rotate(d, dt, about='center'):
            if about == 'center':
                d.rotate(dt*5, about_point=d.get_center())
            elif about == 'start':
                d.rotate(dt*5, about_point=d.get_start())

class RandomOrientations(ThreeDSlide):
    def construct(self):
        setup_all_spins(self)
        self.start_loop()
        self.do_camera_rotation() 
        self.end_loop()

class Spins(ThreeDSlide):
    def construct(self):
        setup_all_spins(self, updater=rotate)
        self.start_loop()
        self.do_camera_rotation() 
        self.end_loop()


class SpinsJoin(ThreeDSlide):
    def construct(self):
        # Start init at the center
        arrows, spheres = setup_all_spins(self, updater=rotate)
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
        self.next_slide()
        
        self.start_loop()
        self.wait(PI)
        self.end_loop()                    
        
        
class SpinsJoinMain(ThreeDSlide):
    def construct(self):
        # Start init at the center
        arrows, spheres = setup_all_spins(self, updater=rotate, locations=np.zeros_like(grid))
        self.remove(*spheres.flatten())
        self.set_camera_orientation(phi=80*DEGREES, theta=45 * DEGREES, zoom=3)
        M0 = get_main_spin()
        self.add(M0)
        self.start_loop()
        do_camera_rotation(self)
        self.end_loop()


class SpinRFPulse(ThreeDSlide):
    def construct(self):
        add_axes(self) 
        self.set_camera_orientation(phi=80*DEGREES, theta=45 * DEGREES, zoom=3)
        M0 = get_main_spin()
        self.add(M0)
        
        self.begin_ambient_camera_rotation(rate=DEGREES*5, about='theta')
        self.play(
            Rotate(M0, PI/2, axis=[1, 1, 0], about_point=M0.get_start()),
            rate_func=exponential_decay,
            run_time=2,
        )
        self.stop_ambient_camera_rotation()
        self.wait(PI)
        
        M0.add_updater(lambda d, dt: rotate(d, dt, about='start'))
        do_camera_rotation(self) 
        