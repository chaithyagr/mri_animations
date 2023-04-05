from manim import *
#from manim.opengl import *
import numpy as np
from manim_slides import ThreeDSlide, Slide

# Configurtaions
config.background_color = WHITE
config.frame_width = 9
config.frame_height = 9
config.pixel_width = 1080
config.pixel_height = 1080

N = 3
grid = np.asarray(np.meshgrid(np.arange(-N//2+1, N//2+1), np.arange(-N//2+1, N//2+1), np.arange(-N//2+1, N//2+1)))
RES = 8
np.random.seed(0)
TEST = False

def get_spin(location, angle_in_xy, angle_in_z=PI/3, radius=0.2, color=RED, sphere_color=BLUE):
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


def get_main_spin(get_trace=False, location=[0, 0, 0], length=1):
    start = location
    end = location + [0, 0, length]
    M0 = Arrow3D(
            start=start,
            end=end,
            resolution=RES,
            color=GREY,
            thickness=0.02,
            height=0.2,
            base_radius=0.05,
        )
    if get_trace:
        for_trace = Arrow(
                start=start,
                end=end + [0, 0, 0.2],
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

def setup_all_spins(obj, updater=None, locations='grid', downwards=False, need_main_spin=False, add_to_scene=True):
    arrows = np.empty((N, N, N), dtype=object)
    spheres =  np.empty((N, N, N), dtype=object)
    traces =  np.empty((N, N, N), dtype=object)
    if locations == 'grid':
        loc = grid
    else:
        loc = locations
    for i in range(N):
        for j in range(N):
            for k in range(N):
                arrow_angle = np.random.random()*2*PI
                arrow_angle_z = PI/3
                if downwards:
                    arrow_angle_z = np.random.random()*2*PI
                if need_main_spin:
                    arrows[i, j, k], spheres[i, j, k] = get_main_spin(
                        location=loc[:, i, j, k],
                        get_trace=True,
                        length=0.6,
                    )
                    traces[i, j, k] = TracedPath(spheres[i, j, k].get_end, stroke_width=4, stroke_color=RED, dissipating_time=PI/4)
                else:
                    arrows[i, j, k], spheres[i, j, k] = get_spin(loc[:, i, j, k], arrow_angle, radius=0.1, angle_in_z=arrow_angle_z)
                if updater is not None:
                    arrows[i, j, k].add_updater(updater)
    if add_to_scene:
        obj.add(*arrows.flatten(), *spheres.flatten())
    if need_main_spin:
        return arrows, spheres, traces
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


def update_relaxation(d, dt, speed=20):
    d.rotate(dt*speed, axis=[0, 0, 1], about_point=d[0].get_start())


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


def do_relax_varying(obj, M0, get_animations=False):
    animations = [
        Rotate(M0, -PI/2, axis=[1, -1, 0], about_point=M0[0].get_start()),
        UpdateFromAlphaFunc(M0, lambda d, dt: update_relaxation(d, dt, speed=(np.random.random()-0.5)*4+15)),
    ]
    if not get_animations:
        obj.play(*animations, run_time=PI)
    else:
        return animations   
     
   
def do_relax_varying(obj, M0, get_animations=False):
    animations = [
        Rotate(M0, -PI/2, axis=[1, -1, 0], about_point=M0[0].get_start()),
        UpdateFromAlphaFunc(M0, lambda d, dt: update_relaxation(d, dt, speed=(np.random.random()-0.5)*4+15)),
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
        
        trace = TracedPath(Trace_M0.get_end, stroke_width=4, stroke_color=RED, dissipating_time=PI/4, n_points_per_cubic_curve=2)
        self.add(trace)
        flip_rf(self, M0)
        do_relax(self, M0)
        
        
        flip_animations = flip_rf(self, M0, get_animations=True)
        relax_animations = do_relax(self, M0, get_animations=True)
        
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
            point[2] = 0
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
            stroke_color=BLUE,
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



class FID3DSplit(ThreeDSlide):
    def construct(self):
        add_axes(self) 
        M0, Trace_Arrow, Trace_M0 = setup_all_spins(self, locations=np.zeros_like(grid), need_main_spin=True, add_to_scene=False)
        self.add(*M0.flatten())
        self.set_camera_orientation(phi=0*DEGREES, theta=90 * DEGREES, zoom=1)
        
        Animations1 = []
        Animations2 = []
        flip_rf_animations = []
        flip_relax_animations = []
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    arrow = M0[i, j, k]
                    end_loc = grid[:, i,j,k]
                    mid_loc = end_loc.copy()
                    mid_loc[2] = 0.5
                    line1 = Line(start=[0,0,0], end=mid_loc)
                    line2 = Line(start=mid_loc, end=end_loc)
                    Animations1.append(MoveAlongPath(arrow, line1))
                    Animations2.append(MoveAlongPath(arrow, line2))
                    flip_rf_animations += flip_rf(self, M0[i,j,k], get_animations=True)
                    flip_relax_animations += do_relax(self, M0[i,j,k], get_animations=True)
                    
        self.play(
            *Animations1, 
            self.camera.animate.set_theta(40*DEGREES),
            run_time=1,
        )
        self.play(
            *Animations2, 
            self.camera.animate.set_phi(140*DEGREES),
            run_time=PI,
        )
        self.move_camera(zoom=0.8, run_time=0.5)


from functools import partial


def get_fid_relax_function(time, get_end, start, inphase=True):
    point = get_end()
    if inphase:
        point[1] = start[1] + (4 + (time-1)/5)
    else:
        point[0] = start[0] + (4 + (time-1)/5)
    point[2] = start[2]
    return point

class FID3D(ThreeDSlide):
    def construct(self):
        add_axes(self) 
        M0, Trace_Arrow, Trace_M0 = setup_all_spins(self, locations=grid, need_main_spin=True, add_to_scene=False)
        self.add(*M0.flatten(), *Trace_M0.flatten())
        self.set_camera_orientation(phi=60*DEGREES, theta=140*DEGREES, zoom=0.8)
        
        flip_rf_animations = []
        flip_relax_animations = []
        inphase_traces = []
        outphase_traces = []
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    flip_rf_animations += flip_rf(self, M0[i,j,k], get_animations=True)
                    flip_relax_animations += do_relax(self, M0[i,j,k], get_animations=True)
                    arrow = M0[i, j, k][0]
                    inphase_traces.append(TracedPath(
                        partial(get_fid_relax_function, get_end=Trace_Arrow[i,j,k].get_end, start=arrow.get_start()),
                        stroke_width=3,
                        stroke_color=BLUE,
                        dissipating_time=PI/4,
                        update_time=True
                    ))
                    outphase_traces.append(TracedPath(
                        partial(get_fid_relax_function, get_end=Trace_Arrow[i,j,k].get_end, start=arrow.get_start(), inphase=False),
                        stroke_width=3,
                        stroke_color=GREEN,
                        dissipating_time=PI/4,
                        update_time=True
                    ))
        self.start_loop()           
        self.play(
            *flip_rf_animations,
            run_time=1,
            rate_function=linear,
        )
        self.add(*inphase_traces, *outphase_traces)
        self.play(
            *flip_relax_animations,
            run_time=PI,
        )
        self.end_loop()

def do_relax_varying(obj, M0, get_animations=False):
    animations = [
        Rotate(M0, -PI/2, axis=[1, -1, 0], about_point=M0[0].get_start()),
        UpdateFromAlphaFunc(M0, lambda d, dt: update_relaxation(d, dt, speed=np.random.uniform(10, 30))),
    ]
    if not get_animations:
        obj.play(*animations, run_time=PI)
    else:
        return animations   

class FID3DGrads(ThreeDSlide):
    def construct(self):
        add_axes(self) 
        M0, Trace_Arrow, Trace_M0 = setup_all_spins(self, locations=grid, need_main_spin=True, add_to_scene=False)
        self.add(*M0.flatten(), *Trace_M0.flatten())
        self.set_camera_orientation(phi=60*DEGREES, theta=140*DEGREES, zoom=0.8)
        
        flip_rf_animations = []
        flip_relax_animations = []
        inphase_traces = []
        outphase_traces = []
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    flip_rf_animations += flip_rf(self, M0[i,j,k], get_animations=True)
                    flip_relax_animations += do_relax_varying(self, M0[i,j,k], get_animations=True)
                    arrow = M0[i, j, k][0]
                    inphase_traces.append(TracedPath(
                        partial(get_fid_relax_function, get_end=Trace_Arrow[i,j,k].get_end, start=arrow.get_start()),
                        stroke_width=3,
                        stroke_color=BLUE,
                        dissipating_time=PI/4,
                        update_time=True
                    ))
                    outphase_traces.append(TracedPath(
                        partial(get_fid_relax_function, get_end=Trace_Arrow[i,j,k].get_end, start=arrow.get_start(), inphase=False),
                        stroke_width=3,
                        stroke_color=GREEN,
                        dissipating_time=PI/4,
                        update_time=True
                    ))
                    
        self.start_loop()           
        self.play(
            *flip_rf_animations,
            run_time=1,
            rate_function=linear,
        )
        self.add(*inphase_traces, *outphase_traces)
        self.play(
            *flip_relax_animations,
            run_time=PI,
        )
        self.end_loop()

def join_acquired_signal(time, Trace_Arrow, inphase=True):
    if inphase:
        point = [0, 3.5, 1]
    else:
        point = [3.5, 0, 1]
    for i in range(N):
        for j in range(N):
            for k in range(N):
                arrow = Trace_Arrow[i, j, k]
                if inphase:
                    point[0] += arrow.get_end()[0]
                    point[1] += time/2
                else:
                    point[1] += arrow.get_end()[1]
                    point[0] += time/2
    return point 
    

class FID3DGradsJoinSig(ThreeDSlide):
    def construct(self):
        add_axes(self) 
        M0, Trace_Arrow, Trace_M0 = setup_all_spins(self, locations=grid, need_main_spin=True, add_to_scene=False)
        self.add(*M0.flatten(), *Trace_M0.flatten())
        self.set_camera_orientation(phi=60*DEGREES, theta=140*DEGREES, zoom=0.8)
        
        flip_rf_animations = []
        flip_relax_animations = []
        inphase_traces = []
        outphase_traces = []
        join_animations = []
        fade_animations = []
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    flip_rf_animations += flip_rf(self, M0[i,j,k], get_animations=True)
                    flip_relax_animations += do_relax_varying(self, M0[i,j,k], get_animations=True)
                    arrow = M0[i, j, k][0]
                    inphase_traces.append(TracedPath(
                        partial(join_acquired_signal, get_end=Trace_Arrow[i,j,k].get_end, start=arrow.get_start()),
                        stroke_width=3,
                        stroke_color=BLUE,
                        dissipating_time=PI/4,
                        update_time=True
                    ))
                    outphase_traces.append(TracedPath(
                        partial(join_acquired_signal, get_end=Trace_Arrow[i,j,k].get_end, start=arrow.get_start(), inphase=False),
                        stroke_width=3,
                        stroke_color=GREEN,
                        dissipating_time=PI/4,
                        update_time=True
                    ))
                    start_loc_inphase = get_fid_relax_function(0, get_end=Trace_Arrow[i,j,k].get_end, start=arrow.get_start())
                    start_loc_outphase = get_fid_relax_function(0, get_end=Trace_Arrow[i,j,k].get_end, start=arrow.get_start(), inphase=False)
                    line_inphase = Line(start=start_loc_inphase, end=[0,3.5,1])
                    line_outphase = Line(start=start_loc_outphase, end=[3.5, 0, 1])
                    join_animations.append(MoveAlongPath(inphase_traces[-1], line_inphase))
                    join_animations.append(MoveAlongPath(outphase_traces[-1], line_outphase))
                    fade_animations.append(FadeOut(inphase_traces[-1]))
                    fade_animations.extend(FadeOut(outphase_traces[-1]))
        inphase_joined = TracedPath(
            partial(get_fid_relax_function, Trace_Arrow=Trace_Arrow),
            stroke_width=5,
            stroke_color=GREEN,
            dissipating_time=PI/4,
            update_time=True
        )
        outphase_joined = TracedPath(
            partial(get_fid_relax_function, Trace_Arrow=Trace_Arrow, inphase=False),
            stroke_width=5,
            stroke_color=GREEN,
            dissipating_time=PI/4,
            update_time=True
        )
        
        self.start_loop()           
        self.play(
            *flip_rf_animations,
            run_time=1,
            rate_function=linear,
        )
        self.add(*inphase_traces, *outphase_traces)
        self.play(
            *flip_relax_animations,
            *join_animations,
            *fade_animations,
            FadeIn(inphase_joined),
            FadeIn(outphase_joined),
            run_time=PI,
        )
        self.end_loop()



