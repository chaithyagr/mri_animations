from manim import *
#from manim.opengl import *
import numpy as np
from manim_slides import ThreeDSlide, Slide

# Configurtaions
config.background_color = WHITE
#config.frame_width = 9
#config.frame_height = 9
#config.pixel_width = 480
#config.pixel_height = 480

N = 3
grid = np.asarray(np.meshgrid(np.arange(-N//2+1, N//2+1), np.arange(-N//2+1, N//2+1), np.arange(-N//2+1, N//2+1)))
RES = 8
np.random.seed(0)
TEST = False

def get_spin(location, angle_in_xy, angle_in_z=PI/3, radius=0.2, color=RED, sphere_color=BLUE, downwards=False):
    sphere = Sphere(
                center=location,
                radius=radius,
                resolution=(RES, RES),
                u_range=[0.001, PI - 0.001],
                v_range=[0, TAU],
            )
    sphere.set_color(sphere_color)
    arrow_radius = radius * 5
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
                thickness=0.002,
                height=0.1,
                base_radius=0.04,
        )
    return arrow, sphere


def get_main_spin(get_trace=False):
    M0 = Arrow3D(
            start=[0, 0, 0],
            end=[0, 0, 1],
            resolution=RES,
            color=BLUE,
            thickness=0.02,
            height=0.2,
            base_radius=0.05,
        )
    if get_trace:
        for_trace = Arrow(
                start=[0, 0, 0],
                end=[0, 0, 1.2],
                resolution=RES,
                color=WHITE,
                thickness=0.002,
                stroke_width=0.001,
            )
        M0 = VGroup(M0, for_trace)
        return M0, for_trace
    else:
        return M0

def do_camera_rotation(obj):
    if TEST == False:
        obj.begin_ambient_camera_rotation(rate=DEGREES*10, about='theta')
        obj.wait(DEGREES*90)
        obj.stop_ambient_camera_rotation()
        obj.begin_ambient_camera_rotation(rate=-DEGREES*10, about='phi')
        obj.wait(DEGREES*90)
        obj.stop_ambient_camera_rotation()
        
        obj.begin_ambient_camera_rotation(rate=-DEGREES*10, about='theta')
        obj.wait(DEGREES*90)
        obj.stop_ambient_camera_rotation()
        obj.begin_ambient_camera_rotation(rate=DEGREES*10, about='phi')
        obj.wait(DEGREES*90)
        obj.stop_ambient_camera_rotation()
    else:
        obj.wait(PI)


def add_axes(obj):
    axes = ThreeDAxes()
    axes.set_color(BLACK)
    obj.add(axes)
    obj.renderer.camera.light_source.move_to(3*IN) # changes the source of the light
    obj.set_camera_orientation(phi=80*DEGREES, theta=45 * DEGREES)

def setup_all_spins(obj, updater=None, locations='grid', downwards=False):
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
                if downwards:
                    arrow_angle_z = np.random.random()*PI
                arrows[i, j, k], spheres[i, j, k] = get_spin(loc[:, i, j, k], arrow_angle, radius=0.1, angle_in_z=arrow_angle_z)
                if updater is not None:
                    arrows[i, j, k].add_updater(updater)
    obj.add(*arrows.flatten(), *spheres.flatten())
    return arrows, spheres

def rotate(d, dt, about='center'):
        if about == 'center':
            d.rotate(dt*5, about_point=d.get_center())
        elif about == 'start':
            d.rotate(dt*5, axis=[0, 0, 1], about_point=d.get_start())
 
class RandomOrientations(ThreeDSlide):
    def construct(self):
        add_axes(self)
        arrows, spheres = setup_all_spins(self, downwards=True)
        self.start_loop()
        do_camera_rotation(self)
        self.end_loop()

class Spins(ThreeDSlide):
    def construct(self):
        add_axes(self)
        arrows, spheres = setup_all_spins(self, updater=rotate)
        self.start_loop()
        do_camera_rotation(self)
        self.end_loop()


class SpinsJoin(ThreeDSlide):
    def construct(self):
        # Start init at the center
        add_axes(self)
        arrows, spheres = setup_all_spins(self, updater=rotate)
        self.play(*[FadeOut(sphere) for sphere in spheres.flatten()])
        self.set_camera_orientation(phi=80*DEGREES, theta=45 * DEGREES, zoom=1)
        
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
        add_axes(self)
        arrows, spheres = setup_all_spins(self, updater=rotate, locations=np.zeros_like(grid))
        self.remove(*spheres.flatten())
        self.set_camera_orientation(phi=80*DEGREES, theta=45 * DEGREES, zoom=3)
        M0 = get_main_spin()
        self.add(M0)
        self.start_loop()
        do_camera_rotation(self)
        self.end_loop()


def update_relaxation(d, dt):
    d.rotate(dt*20, axis=[0, 0, 1], about_point=d[0].get_start())


def flip_rf(obj, M0, get_animations=False):
    animations = [
        Rotate(M0, PI/2, axis=[1, -1, 0], about_point=M0[0].get_start()),
    ]
    if not get_animations:
        obj.play(*animations, run_time=1, rate_func=exponential_decay,)
    else:
        return animations
            
def do_relax(obj, M0, get_animations=False):
    animations = [
        Rotate(M0, -PI/2, axis=[1, -1, 0], about_point=M0[0].get_start()),
        UpdateFromAlphaFunc(M0, update_relaxation),
    ]
    if not get_animations:
        obj.play(*animations, run_time=PI)
    else:
        return animations   
    
class SpinRFPulse(ThreeDSlide):
    def construct(self):
        add_axes(self) 
        
        self.set_camera_orientation(phi=80*DEGREES, theta=45 * DEGREES, zoom=3)
        M0, Trace_M0 = get_main_spin(get_trace=True)
        self.add(M0)
        
        flip_rf(self, M0)
        do_relax(self, M0)
        
        trace = TracedPath(Trace_M0.get_end, stroke_width=4, stroke_color=RED, dissipating_time=PI/4, n_points_per_cubic_curve=2)
        self.add(trace)
        flip_rf(self, M0)
        do_relax(self, M0)
        
        
        flip_animations = flip_rf(self, M0, get_animations=True)
        relax_animations = do_relax(self, M0, get_animations=True)
        
        self.play(
            *flip_animations,
            self.camera.animate.set_theta(60*DEGREES),
            run_time=1,
        )
        self.play(
            *relax_animations,
            self.camera.animate.set_phi(60*DEGREES),
            run_time=PI,
        )
        
        self.play(
            *flip_animations,
            self.camera.animate.set_theta(90*DEGREES),
            run_time=1,
        )
        self.play(
            *relax_animations,
            self.camera.animate.set_phi(0*DEGREES),
            run_time=PI,
        )
        self.play(self.camera.animate.set_theta(90*DEGREES), run_time=1)
        self.next_slide() 
        
        self.start_loop()
        flip_rf(self, M0)
        do_relax(self, M0)
        self.end_loop()


class SpinRFPulseCoil(ThreeDSlide):
    def construct(self):
        def get_relax_function(get_end, time, inphase=True):
            point = get_end()
            if inphase:
                point[1] = 1 + (time-1)/2
            else:
                point[0] = 1 + (time-1)/2
            return point
            
            
        add_axes(self) 
        M0, Trace_M0 = get_main_spin(get_trace=True)
        self.add(M0)
        self.set_camera_orientation(phi=0*DEGREES, theta=90 * DEGREES, zoom=3)
        
        trace = TracedPath(Trace_M0.get_end, stroke_width=4, stroke_color=RED, dissipating_time=PI/4)
        self.add(trace)
        
        self.move_camera(phi=0*DEGREES, theta=90 * DEGREES, zoom=1, run_time=1)
        trace_inphase = TracedPath(
            lambda time: get_relax_function(Trace_M0.get_end, time=time),
            stroke_width=5,
            stroke_color=RED,
            dissipating_time=PI/4,
            update_time=True
        )
        trace_outphase = TracedPath(
            lambda time: get_relax_function(Trace_M0.get_end, time=time, inphase=False),
            stroke_width=5,
            stroke_color=GREEN,
            dissipating_time=PI/4,
            update_time=True
        )
        self.next_slide() 
        
        self.start_loop()
        flip_rf(self, M0)
        self.add(trace_inphase, trace_outphase)
        do_relax(self, M0)
        self.end_loop()



class RFPulse(Scene):
    def construct(self):
        # Create a circle to represent the RF pulse
        rf_pulse = Circle(radius=2, fill_opacity=0.5, fill_color=BLUE, stroke_width=0)
        # Add a label for the RF pulse
        rf_label = Text("RF Pulse", font_size=60, color=WHITE).next_to(rf_pulse, DOWN)
        # Create an arrow to represent the magnetization vector
        magnetization_vector = Arrow(start=ORIGIN, end=[0, 2, 0], color=RED, buff=0, stroke_width=10)
        # Add a label for the magnetization vector
        magnetization_label = Text("Magnetization", font_size=40, color=RED).next_to(magnetization_vector, RIGHT)
        
        # Add the RF pulse and magnetization vector to the scene and animate them
        self.play(Create(rf_pulse), Create(magnetization_vector), run_time=2*PI)
        # Apply a rotation animation to the magnetization vector
        self.play(Rotate(magnetization_vector, angle=TAU/4, about_point=ORIGIN, run_time=0.5))
        # Add the labels to the scene and animate them
        self.play(Create(rf_label), Create(magnetization_label))
        # Remove the RF pulse, magnetization vector, and labels from the scene
        self.play(FadeOut(rf_pulse), FadeOut(magnetization_vector), FadeOut(rf_label), FadeOut(magnetization_label))
