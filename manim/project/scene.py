import math
import random
import numpy as np
from manim import *

CITATION_SIZE = 0.3


class Intro(Scene):
    def construct(self):
        # Greetings and Justin's video screenshot
        hi = Text("Hey everyone, I'm James.")
        justin = ImageMobject("images/Justin's Video.png")
        hi.to_edge(UP, buff=0.2)
        justin.scale(0.9)
        justin.next_to(hi, DOWN)
        self.play(Write(hi))
        self.wait(9)

        self.play(FadeIn(justin))
        self.wait(15)

        self.play(FadeOut(hi), FadeOut(justin))

        # Problem to discuss in this video
        problem = Text(
            "How to find a uniform random point\nin a high-dimensional ball?",
            t2c={"high-dimensional ball": RED})
        problem.shift(UP)
        circ_def = MathTex(
            r"C(r) = \{(x, y) : x^2 + y^2 = r^2\}")
        ball_def = MathTex(
            r"B_n(r) = \{(x_1, \ldots, x_n) : x_1^2 + \ldots + x_n^2 = r^2\}")
        circ_def.next_to(problem, DOWN)
        ball_def.next_to(circ_def, DOWN)

        self.play(Create(problem), run_time=2)
        self.wait(15)

        self.play(Create(circ_def))
        self.wait(7)

        self.play(Create(ball_def))
        self.wait(8)

        # Fade out all Mobjects
        self.play(*[FadeOut(obj) for obj in self.mobjects])


class RejSampling(Scene):
    def construct(self):
        # Title
        rej_samp_title = Text("Rejection Sampling")
        rej_samp_title.to_edge(UP, buff=0.2)
        self.play(Create(rej_samp_title))
        self.wait(5)

        # Rejection sampling animation
        r = 3
        square = Square(side_length=r * 2)
        circle = Circle(radius=r)

        self.play(Create(square), Create(circle))
        self.wait(1)

        num_points = 10  # Must be greater than 1
        assert num_points > 1
        points = [[0.39 * r, -0.57 * r, 0], [-0.73 * r, 0.82 * r, 0]]

        for _ in range(num_points - 2):
            points.append([2 * (random.random() - 0.5) * r,
                           2 * (random.random() - 0.5) * r,
                           0])

        for p in points:
            dot = Dot(point=p)
            self.add(dot)

            if (p[0] ** 2 + p[1] ** 2 <= r ** 2):
                self.play(dot.animate.set_color(GREEN))
            else:
                self.play(dot.animate.set_color(RED))
            self.wait(0.5)

        self.wait(1)
        self.play(*[FadeOut(obj) for obj in self.mobjects])
        self.wait(1)

        # Code
        rs_code_raw = '''
        def rejection_sampling(n):
            trials = 0
            while True:
                trials += 1
                vec = np.empty(n)
                for i in range(n):
                    vec[i] = (random.random() - 0.5) * 2
                if (np.dot(vec, vec) <= 1):
                    return (vec, trials)
        '''
        rs_code_rendered = Code(
            code=rs_code_raw,
            language="Python",
            font="Monospace",
            background_stroke_width=0,
            insert_line_no=False)
        self.play(FadeIn(rs_code_rendered))
        self.wait(35)

        # Results
        trials_linear = ImageMobject(
            "images/Trials vs. Dimension Linear.png")
        trials_linear.scale(0.8)
        trials_linear.align_to(ORIGIN - [0.2, 0, 0], RIGHT)

        trials_log = ImageMobject(
            "images/Trials vs. Dimension Log.png")
        trials_log.scale(0.8)
        trials_log.align_to(ORIGIN + [0.2, 0, 0], LEFT)

        self.play(FadeOut(rs_code_rendered))
        self.play(FadeIn(trials_linear))
        self.wait(28)

        self.play(FadeIn(trials_log))
        self.wait(12)

        self.play(*[FadeOut(obj) for obj in self.mobjects])


class BallVolumeIntro(Scene):
    def construct(self):
        x_trials = Text(
            "X: Number of trials to select a valid point", t2s={"X": ITALIC})
        x_trials.scale(0.7)
        x_geom = MathTex(r"X \sim \mathrm{Geom}(p)")
        beta_n = MathTex(r"\beta_n = \mathrm{Vol(Ball)}")
        kappa_n = MathTex(r"\kappa_n = \mathrm{Vol(Box)} = 2^n")
        prob = MathTex(
            r"p = P(\text{Random point in ball}) = \frac{\mathrm{Vol(Ball)}}{\mathrm{Vol(Box)}} = \frac{\beta_n}{\kappa_n}")
        txt_group = VGroup(x_trials, x_geom, beta_n, kappa_n, prob)
        txt_group.arrange(DOWN, aligned_edge=LEFT)

        self.play(FadeIn(x_trials))
        self.wait(9)

        self.play(FadeIn(x_geom))
        self.wait(3)

        self.play(FadeIn(beta_n))
        self.wait(2.5)

        self.play(FadeIn(kappa_n))
        self.wait(19)

        self.play(FadeIn(prob))
        self.wait(9)

        self.play(Indicate(beta_n))
        self.wait(2)

        self.play(*[FadeOut(obj) for obj in self.mobjects])


class BallIntegral(ThreeDScene):
    def construct(self):
        citation = Text("Adapted from Vector Calculus, Linear Algebra, and Differential Forms: A Unified Approach by John Hubbard and Barbara Burke Hubbard", t2s={
            "Vector Calculus, Linear Algebra, and Differential Forms: A Unified Approach": ITALIC})
        citation.scale(CITATION_SIZE)
        citation.to_edge(DOWN + RIGHT)
        self.add_fixed_in_frame_mobjects(citation)

        axes = ThreeDAxes(x_range=(-3.6, 3.6, 0.6),
                          y_range=(-3.6, 3.6, 0.6),
                          z_range=(-2.4, 2.4, 0.6),
                          x_length=12,
                          y_length=12,
                          z_length=8)
        self.set_camera_orientation(
            phi=80 * DEGREES, theta=-95 * DEGREES, distance=5)
        self.add(axes)

        r = 3

        # sphere = Surface(lambda u, v: np.array([
        #     r * np.sin(u) * np.cos(v),
        #     r * np.sin(u) * np.sin(v),
        #     r * np.cos(u)
        # ]), u_range=(0, PI), v_range=(0, 2 * PI))

        # Create the sphere
        sphere = Sphere(ORIGIN, radius=r)
        sphere.set_opacity(0.85)
        self.play(Create(sphere))
        self.wait(3)

        # Dim the sphere
        self.play(sphere.animate.set_opacity(0.15))

        def circle_fun(x):
            ''' Generates the parametric function of a circle with the given
                x coordinate on the sphere with radius r.
            '''
            return ParametricFunction(lambda t: np.array([
                x,
                np.sqrt(r ** 2 - x ** 2) * np.cos(t),
                np.sqrt(r ** 2 - x ** 2) * np.sin(t)
            ]), t_range=(0, 2 * PI))

        def strip_fun(x, dx):
            ''' Generates the parametric function of a thin strip surface
                between the circle at x and the circle at x + dx.
            '''
            return Surface(lambda u, v: np.array([
                x + u * dx,
                np.sqrt(r ** 2 - (x + u * dx) ** 2) * np.cos(v),
                np.sqrt(r ** 2 - (x + u * dx) ** 2) * np.sin(v)
            ]), u_range=(0, 1), v_range=(0, 2 * PI))

        x_0 = r * 0.4
        dx = r * 0.1
        x_1 = x_0 + dx

        # Create the circles and the thin strip
        x_0_circ = circle_fun(x_0)
        x_1_circ = circle_fun(x_1)
        x_0_circ.set_color(YELLOW)
        x_1_circ.set_color(GREEN)

        strip = strip_fun(x_0, dx)
        strip.set_color(RED)
        self.play(Create(x_0_circ), Create(x_1_circ), Create(strip))
        self.wait(26)

        # Draw the right triangle and label each line segment
        r_line = Line3D(start=ORIGIN, end=np.array(
            [x_0, 0, np.sqrt(r ** 2 - x_0 ** 2)]), thickness=0.02, color=WHITE)
        r_label = MathTex(r"1")
        r_label.rotate(PI / 2, axis=RIGHT)
        r_label.move_to(r_line.get_center() + UL * MED_SMALL_BUFF)

        x_line = Line3D(start=ORIGIN, end=np.array(
            [x_0, 0, 0]), thickness=0.02, color=WHITE)
        x_label = MathTex(r"x_1")
        x_label.rotate(PI / 2, axis=RIGHT)
        x_label.next_to(x_line, UP, buff=LARGE_BUFF)

        circ_r_line = Line3D(start=np.array([x_0, 0, 0]), end=np.array(
            [x_0, 0, np.sqrt(r ** 2 - x_0 ** 2)]), thickness=0.02, color=WHITE)
        circ_r_label = MathTex(r"r")
        circ_r_label.rotate(PI / 2, axis=RIGHT)
        circ_r_label.next_to(circ_r_line, RIGHT, buff=SMALL_BUFF)

        self.play(Create(r_line))
        self.play(Write(r_label))
        self.wait(1)

        self.play(Create(x_line))
        self.play(Write(x_label))
        self.wait(1)

        self.play(Create(circ_r_line))
        self.play(Write(circ_r_label))
        self.wait(1)

        # Volume of thin strip
        r_expr = MathTex(r"r = \sqrt{1 - x_1^2}")
        dV = MathTex(r"dV = r^{n-1} \beta_{n-1} \, dx_1")
        txt_group = VGroup(r_expr, dV)
        txt_group.arrange(DOWN, aligned_edge=LEFT)
        txt_group.to_edge(UP + RIGHT)
        self.add_fixed_in_frame_mobjects(txt_group)
        self.remove(txt_group)
        self.play(Create(txt_group[0]))
        self.wait(20)
        self.play(Create(txt_group[1]))
        self.wait(7)

        self.play(*[FadeOut(obj) for obj in self.mobjects])


class BallVolume(Scene):
    def construct(self):
        citation = Text("Adapted from Vector Calculus, Linear Algebra, and Differential Forms: A Unified Approach by John Hubbard and Barbara Burke Hubbard", t2s={
            "Vector Calculus, Linear Algebra, and Differential Forms: A Unified Approach": ITALIC})
        citation.scale(CITATION_SIZE)
        citation.to_edge(DOWN + RIGHT)
        self.add(citation)

        # Relation between beta_n and beta_{n-1}
        r_x1 = MathTex(r"r = \sqrt{1 - x_1^2}")
        beta_n = MathTex(
            r"\beta_n &= \int dV = \int_{-1}^1 r^{n-1} \beta_{n-1} \, dx_1")
        beta_n_2 = MathTex(
            r"\beta_n &= \beta_{n-1} \int_{-1}^1 (1 - x_1^2)^{\frac{n-1}{2}} \, dx_1")
        r_x1.next_to(beta_n, UP)

        self.play(FadeIn(r_x1))
        self.play(Create(beta_n), run_time=2)
        self.wait(6)

        self.play(ReplacementTransform(beta_n, beta_n_2), FadeOut(r_x1))
        self.wait(3)

        # Multiplying factor c_n
        t_int = MathTex(r"c_n = \int_{-1}^1 (1 - t^2)^{\frac{n-1}{2}} \, dt")
        t_int.save_state()
        t_int.next_to(beta_n_2, DOWN)
        t_int.scale(0.7)

        beta_cn = MathTex(r"\beta_n = c_n \beta_{n - 1}")
        beta_cn.save_state()
        beta_cn.next_to(t_int, DOWN)
        beta_cn.scale(0.7)

        # It remains to figure out c_n
        self.play(FadeIn(t_int))
        self.wait(5)

        # Recurrence relation
        self.play(FadeIn(beta_cn))
        self.wait(2)

        self.play(FadeOut(beta_n_2),
                  FadeOut(beta_cn),
                  Restore(t_int))

        # Trig substitution
        trig_sub = MathTex(r"t = \sin(\theta)")
        trig_sub.next_to(t_int, DOWN)
        trig_sub.scale(0.7)
        dt = MathTex(r"dt = \cos(\theta) \, d\theta")
        dt.next_to(trig_sub, DOWN)
        dt.scale(0.7)
        self.play(FadeIn(trig_sub))
        self.play(FadeIn(dt))

        trig_int_n = MathTex(
            r"c_n &= \int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} \cos^{n}(\theta) \, d\theta")
        trig_int_n_2 = MathTex(
            r"c_{n-2} &= \int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} \cos^{n-2}(\theta) \, d\theta")
        # trig_ints = MathTex(r"c_n & = \int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} \cos^{n}(\theta) \, d\theta \\ c_{n-2} & = \int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} \cos^{n-2}(\theta) \, d\theta")

        self.play(FadeOut(trig_sub), FadeOut(dt),
                  ReplacementTransform(t_int, trig_int_n))
        self.play(trig_int_n.animate.shift(UP))
        trig_int_n_2.next_to(trig_int_n, DOWN)
        self.play(Create(trig_int_n_2))

        self.play(FadeOut(trig_int_n_2))
        self.play(trig_int_n.animate.shift(DOWN))

        # Integration by parts
        u_dv = MathTex(
            r"u &= \cos^{n-1}(\theta) \\ dv &= \cos(\theta) \, d\theta")
        du_v = MathTex(
            r"du &= -\sin(\theta) \cos^{n-2}(\theta) \, d\theta \\ v &= \sin(\theta)")
        u_dv.scale(0.7)
        du_v.scale(0.7)
        u_dv.next_to(trig_int_n, DOWN)
        u_dv.shift(1.5 * LEFT)
        du_v.next_to(u_dv, 2 * RIGHT)
        self.play(FadeIn(u_dv))
        self.play(FadeIn(du_v))

        ibp = MathTex(
            r"c_n = \cos^{n-1}(\theta) \sin(\theta) \bigg|_{-\frac{\pi}{2}}^{\frac{\pi}{2}} + \int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} \sin^2(\theta) (n-1) \cos^{n-2}(\theta) \, d\theta")
        ibp_text = Text("Integration by Parts", font_size=36)
        ibp_text.next_to(ibp, UP)
        self.play(FadeOut(u_dv), FadeOut(du_v),
                  ReplacementTransform(trig_int_n, ibp), FadeIn(ibp_text))

        # First term on RHS is 0
        cos_is_zero = MathTex(
            r"\cos \left(-\frac{\pi}{2} \right) = \cos{\frac{\pi}{2}} = 0")
        cos_is_zero.next_to(ibp, DOWN)
        cos_is_zero.scale(0.7)
        self.play(FadeOut(ibp_text))
        self.play(FadeIn(cos_is_zero))
        self.wait(3)

        ibp_one_term = MathTex(
            r"c_n = (n-1) \int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} \sin^2(\theta) \cos^{n-2}(\theta) \, d\theta")
        self.play(ReplacementTransform(
            ibp, ibp_one_term), FadeOut(cos_is_zero))

        # Trig identities
        sin_cos_id = MathTex(r"\sin^2(\theta) = 1 - \cos^2(\theta)")
        sin_cos_id.next_to(ibp_one_term, DOWN)
        sin_cos_id.scale(0.7)
        self.play(FadeIn(sin_cos_id))

        ibp_split = MathTex(
            r"c_n = (n-1) \left[\int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} \cos^{n-2}(\theta) \, d\theta - \int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} \cos^n(\theta) \, d\theta\right]")
        self.play(ReplacementTransform(
            ibp_one_term, ibp_split), FadeOut(sin_cos_id))

        self.play(ibp_split.animate.shift(UP))
        trig_int_n = MathTex(
            r"c_n &= \int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} \cos^{n}(\theta) \, d\theta")
        trig_int_n.next_to(ibp_split, DOWN)
        trig_int_n_2.next_to(trig_int_n, DOWN)
        self.play(FadeIn(trig_int_n))
        self.play(FadeIn(trig_int_n_2))

        # Rearranging
        ibp_rearrange = MathTex(r"c_n = (n-1) c_{n-2} - (n-1) c_n")
        self.play(ReplacementTransform(
            ibp_split, ibp_rearrange), FadeOut(trig_int_n), FadeOut(trig_int_n_2))

        # Recurrence relation between c_n and c_{n-2}
        result = MathTex(r"c_n = \frac{n-1}{n}c_{n-2}")
        self.play(ReplacementTransform(ibp_rearrange, result))
        self.wait(5)

        # End of first part of recording

        # Base cases
        base_text = Text("Base cases:")
        c_1 = MathTex(r"c_1 = \int_{-1}^1 (1 - t^2)^0 \, dt = 2")
        c_1.next_to(base_text, DOWN)
        c_2 = MathTex(
            r"c_2 = \int_{-1}^1 (1 - t^2)^{\frac{1}{2}} \, dt = \frac{\pi}{2}")
        c_2.next_to(c_1, DOWN)

        base_group = VGroup(base_text, c_1, c_2)
        base_group.move_to(ORIGIN)
        self.play(result.animate.to_edge(UP))
        self.play(FadeIn(base_group))
        self.wait(6)

        result.generate_target()
        result.target.scale(0.7)
        result.target.to_edge(UP + RIGHT)

        base_short = MathTex(r"c_1 &= 2 \\ c_2 &= \frac{\pi}{2}")
        base_short.next_to(result.target, DOWN)
        base_short.scale(0.7)

        beta_cn.restore()
        beta_cn.shift(2 * UP)
        self.play(ReplacementTransform(base_group, base_short),
                  MoveToTarget(result),
                  FadeIn(beta_cn))
        self.wait(5)

        # Ratio between beta_n and kappa_n
        kappa_rel = MathTex(
            r"\kappa_{n - 1} \stackrel{\times 2}{\longrightarrow} \kappa_n")
        kappa_rel.next_to(beta_cn, DOWN)
        self.play(FadeIn(kappa_rel), FadeOut(citation))
        self.wait(5)

        beta_rel = MathTex(
            r"\beta_{n - 1} \stackrel{\times c_n}{\longrightarrow} \beta_n")
        beta_rel.next_to(kappa_rel, DOWN, buff=MED_LARGE_BUFF)
        self.play(FadeIn(beta_rel))
        self.wait(4)

        bk_group = VGroup(kappa_rel, beta_rel)

        ratio_rel = MathTex(
            r"\frac{\beta_{n - 1}}{\kappa_{n-1}} \stackrel{\times \frac{c_n}{2}}{\longrightarrow} \frac{\beta_n}{\kappa_n}"
        )
        ratio_rel.next_to(beta_cn, DOWN, buff=MED_LARGE_BUFF)
        self.play(FadeOut(bk_group),
                  FadeIn(ratio_rel))
        self.wait(8)

        # End of second part of recording

        self.wait(15)

        prob = MathTex(
            r"P(\text{Random point in ball}) = \frac{\mathrm{Vol(Ball)}}{\mathrm{Vol(Box)}} = \frac{\beta_n}{\kappa_n}")
        prob.next_to(ratio_rel, DOWN)
        self.play(FadeIn(prob))
        self.wait(29)

        ex = MathTex(r"E[X] = \frac{1}{P(\text{Random point in ball)}}")
        incr_fast = Text("increases faster than exponentially", color=RED)
        ex_group = VGroup(ex, incr_fast)
        ex_group.arrange(DOWN)
        ex_group.next_to(prob, DOWN)
        self.play(FadeIn(ex))
        self.wait(11)
        self.play(GrowFromCenter(incr_fast))
        self.wait(1)

        self.play(*[FadeOut(obj) for obj in self.mobjects])


class SphericalCoord(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range=(-3.6, 3.6, 0.6),
                          y_range=(-3.6, 3.6, 0.6),
                          z_range=(-2.4, 2.4, 0.6),
                          x_length=12,
                          y_length=12,
                          z_length=8)
        self.set_camera_orientation(
            phi=70 * DEGREES, theta=60 * DEGREES, distance=5)
        self.add(axes)

        r = 3

        sphere = Surface(lambda u, v: np.array([
            r * np.sin(u) * np.cos(v),
            r * np.sin(u) * np.sin(v),
            r * np.cos(u)
        ]), u_range=(0, PI), v_range=(0, 2 * PI))
        sphere.set_opacity(0.7)

        def circle_fun(phi):
            ''' Generates the parametric function of a circle with the given
                phi coordinate on the sphere with radius r. '''
            return ParametricFunction(lambda t: np.array([
                r * np.sin(phi) * np.cos(t),
                r * np.sin(phi) * np.sin(t),
                r * np.cos(phi)
            ]), t_range=(0, 2 * PI))

        phi_eqt = 90 * DEGREES  # Equator
        phi_top = 20 * DEGREES
        equator = circle_fun(phi_eqt)
        equator.set_color(RED)
        top_circle = circle_fun(phi_top)
        top_circle.set_color(GREEN)

        self.play(Create(sphere))
        self.wait(36)

        self.play(Create(equator))
        self.wait(8)

        self.play(Create(top_circle))
        self.wait(15)

        self.play(*[FadeOut(obj) for obj in self.mobjects])


class SphereDistIntro(Scene):
    def construct(self):
        r = 3
        circle = Circle(radius=r)
        self.play(Create(circle))
        self.wait(2)

        theta_0 = PI / 3
        unit_vec = Vector([r * np.cos(theta_0), r * np.sin(theta_0), 0])
        self.play(Create(unit_vec))
        self.play(Rotate(unit_vec,
                         angle=2 * PI,
                         about_point=ORIGIN,
                         run_time=2,
                         rate_func=smooth))
        self.wait(1)
        self.play(ApplyMethod(unit_vec.put_start_and_end_on,
                              [0, 0, 0],
                              [0.2 * r * np.cos(theta_0),
                               0.2 * r * np.sin(theta_0), 0],
                              run_time=2,
                              rate_func=there_and_back))

        graph_group = VGroup(circle, unit_vec)
        graph_group.generate_target()
        graph_group.target.scale(0.6)
        graph_group.target.to_edge(LEFT)
        self.play(MoveToTarget(graph_group))
        self.wait(1)

        tasks = BulletedList("Choose a uniformly random direction",
                             "Choose a random distance from the origin")
        tasks.scale(0.9)
        rand_dir = Text("Random direction = Random point on sphere")
        rand_dir.scale(0.6)
        sphere_def = MathTex(
            r"S_n(r) = \{(x_1, \ldots, x_n) : x_1^2 + \ldots + x_n^2 = r^2\}")
        txt_group = VGroup(tasks, rand_dir, sphere_def)
        txt_group.arrange(DOWN, buff=MED_LARGE_BUFF)
        txt_group.next_to(graph_group, RIGHT, buff=MED_LARGE_BUFF)

        self.play(FadeIn(tasks[0]))
        self.wait(1)

        self.play(FadeIn(tasks[1]))
        self.wait(1)

        self.play(Indicate(tasks[0]), FadeIn(rand_dir))
        self.wait(8)

        self.play(FadeIn(sphere_def))
        self.wait(20)

        self.play(*[FadeOut(obj) for obj in self.mobjects])


def get_diag_arc(r, eps, color):
    ''' Returns an arc in the upper right corner
        with radius r and width about 2 * r * eps. '''
    diag_arc_start = math.atan2(1 - 2 * eps, 1)
    diag_arc_end = math.atan2(1, 1 - 2 * eps)
    diag_arc = Arc(arc_center=ORIGIN,
                   radius=r,
                   start_angle=diag_arc_start,
                   angle=diag_arc_end - diag_arc_start)
    diag_arc.set_color(color)
    return diag_arc


def get_vert_arc(r, eps, color):
    ''' Returns an arc in the upper edge
        with radius r and width about 2 * r * eps. '''
    vert_arc_start = math.atan2(1, eps)
    vert_arc_end = math.atan2(1, -eps)
    vert_arc = Arc(arc_center=ORIGIN,
                   radius=r,
                   start_angle=vert_arc_start,
                   angle=vert_arc_end - vert_arc_start)
    vert_arc.set_color(color)
    return vert_arc


def play_normalize_vec(scene, arr, r):
    ''' Creates a vector based on arr and then normalizes it.
        Returns the unit vector created.
        Duration: 4 seconds
    '''
    unit_arr = arr / np.linalg.norm(arr) * r
    vec = Vector(arr)
    unit_vec = Vector(unit_arr)
    scene.play(Create(vec))
    scene.wait(1)
    scene.play(ReplacementTransform(vec, unit_vec))
    scene.wait(1)
    return unit_vec  # So that we can fade it out later


class PointOnSphereWrong(Scene):
    def construct(self):
        citation = Text(
            "I learned this result from Professor Robert Kleinberg.")
        citation.scale(CITATION_SIZE)
        citation.to_edge(DOWN + RIGHT)
        self.add(citation)

        r = 3
        square = Square(side_length=r * 2)
        circle = Circle(radius=r)

        self.play(Create(square), Create(circle))
        self.wait(5)

        arrs = [np.array([0.6 * r, -0.4 * r, 0]),
                np.array([-0.93 * r, -0.78 * r, 0])]

        unit_vecs = []

        for arr in arrs:
            unit_vecs.append(play_normalize_vec(self, arr, r))

        self.wait(27)
        self.play(*[FadeOut(unit_vec) for unit_vec in unit_vecs])

        # End of first part of recording

        diag = Line(ORIGIN, [r, r, 0])
        vert = Line(ORIGIN, [0, r, 0])
        diag_label = MathTex(r"\sqrt{2}")
        vert_label = MathTex(r"1")
        diag_label.move_to(diag.get_center() + DR * MED_SMALL_BUFF)
        vert_label.next_to(vert, buff=SMALL_BUFF)
        diag_group = VGroup(diag, diag_label)
        vert_group = VGroup(vert, vert_label)

        self.play(Create(diag_group))
        self.wait(4)
        self.play(Create(vert_group))
        self.wait(4)
        self.play(FadeOut(diag_group), FadeOut(vert_group))

        eps = 0.1
        diag_poly = Polygon([r * (1 - 2 * eps), r, 0],
                            [r, r * (1 - 2 * eps), 0],
                            ORIGIN)
        diag_poly.set_fill(GREEN)
        diag_poly.set_opacity(0.6)
        vert_poly = Polygon([-r * eps, r, 0], [r * eps, r, 0], ORIGIN)
        vert_poly.set_fill(BLUE)
        vert_poly.set_opacity(0.6)

        self.play(FadeIn(diag_poly))
        self.wait(1)
        self.play(FadeIn(vert_poly))
        self.wait(1)

        diag_arc = get_diag_arc(r, eps, GREEN_D)
        diag_arr = np.array([0.87 * r, (0.87 - 0.2 * eps) * r, 0])
        play_normalize_vec(self, diag_arr, r)

        self.play(Create(diag_arc))
        self.play(Indicate(diag_arc))
        self.wait(1)

        vert_arc = get_vert_arc(r, eps, BLUE_D)
        vert_arr = np.array([-eps * 0.05 * r, 0.47 * r, 0])
        play_normalize_vec(self, vert_arr, r)

        self.play(Create(vert_arc))
        self.play(Indicate(vert_arc))
        self.wait(0.5)

        graph_group = Group(*[obj for obj in self.mobjects])
        graph_group.remove(citation)
        self.play(graph_group.animate.to_edge(LEFT))

        len_eq = Text("Length of green arc = Length of blue arc",
                      t2c={"green arc": GREEN_D, "blue arc": BLUE_D})
        prob_ratio = MathTex(
            r"\frac{P(\text{Point falls on green arc})}{P(\text{Point falls on blue arc})}")
        prob_eq_area = MathTex(
            r"= \frac{\mathrm{Area}(\text{Green cone})}{\mathrm{Area}(\text{Blue cone})}")
        prob_eq_vol = MathTex(
            r"= \frac{\mathrm{Vol}(\text{Green cone})}{\mathrm{Vol}(\text{Blue cone})}")
        root_2_sqr = MathTex(r"=\left( \frac{\sqrt{2}}{1} \right)^2")
        root_n_pow = MathTex(r"=\left( \frac{\sqrt{n}}{1} \right)^n")
        incr_fast = Text("increases faster than exponentially", color=RED)

        txt_group = VGroup(len_eq, prob_ratio, prob_eq_area,
                           root_2_sqr, incr_fast)
        txt_group.scale(0.6)
        txt_group.arrange(DOWN, aligned_edge=LEFT)
        txt_group.next_to(graph_group)
        root_n_pow.scale(0.6)
        prob_eq_vol.scale(0.6)
        prob_eq_vol.move_to(prob_eq_area)
        root_n_pow.move_to(root_2_sqr)

        self.play(FadeIn(len_eq))
        self.wait(1)

        self.play(FadeIn(prob_ratio))
        self.wait(1)

        self.play(FadeIn(prob_eq_area))
        self.wait(3)

        self.play(FadeIn(root_2_sqr))
        self.wait(10)

        self.play(ReplacementTransform(root_2_sqr, root_n_pow),
                  ReplacementTransform(prob_eq_area, prob_eq_vol))
        self.wait(5)

        self.play(GrowFromCenter(incr_fast))
        self.wait(15)

        self.play(*[FadeOut(obj) for obj in self.mobjects])


class RotateSquare(Scene):
    def construct(self):
        r = 2
        square = Square(side_length=r * 2)
        self.play(Create(square))
        self.wait(6)
        self.play(Rotate(square, angle=PI / 4))
        self.wait(12)
        self.play(Uncreate(square))


class RotatePlane(LinearTransformationScene):
    def __init__(self):
        LinearTransformationScene.__init__(
            self
        )

    def construct(self):
        q_mat_2d = MathTex(
            r"Q = \begin{bmatrix} \cos(\theta) & -\sin(\theta) \\ \sin(\theta) & \cos(\theta) \end{bmatrix}")

        q_mat_2d.to_corner(UR)
        self.add_foreground_mobject(q_mat_2d)  # Foreground mobjects don't move

        theta = PI / 3
        matrix = [[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta), np.cos(theta)]]
        self.apply_matrix(matrix)
        self.wait(12)

        self.play(*[FadeOut(obj) for obj in self.mobjects])


class PointOnSphereRight(Scene):
    def construct(self):
        # Rotation and orthogonal matrices
        q_mat_2d = MathTex(
            r"Q = \begin{bmatrix} \cos(\theta) & -\sin(\theta) \\ \sin(\theta) & \cos(\theta) \end{bmatrix}")

        orthonormal = MathTex(
            r"\begin{bmatrix} \cos(\theta) \\ \sin(\theta) \end{bmatrix} \text{ and } \begin{bmatrix} -\sin(\theta) \\ \cos(\theta) \end{bmatrix} \text{ are orthonormal}")
        norm_1 = MathTex(r"\sqrt{\cos^2(\theta) + \sin^2(\theta)} = 1")
        norm_1.scale(0.7)
        dot_0 = MathTex(
            r"\begin{bmatrix} \cos(\theta) \\ \sin(\theta) \end{bmatrix} \cdot \begin{bmatrix} -\sin(\theta) \\ \cos(\theta) \end{bmatrix} = 0")
        dot_0.scale(0.7)

        rows_also = Text(
            "Similarly, the row vectors of Q are also orthonormal.", t2s={"Q": ITALIC})
        rows_also.scale(0.6)

        q_group = VGroup(q_mat_2d, orthonormal, norm_1, dot_0, rows_also)
        q_group.arrange(DOWN)

        self.play(FadeIn(q_mat_2d))
        self.wait(5)

        self.play(FadeIn(orthonormal))
        self.wait(4)

        self.play(FadeIn(norm_1))
        self.wait(2)

        self.play(FadeIn(dot_0))
        self.wait(3)

        self.play(FadeIn(rows_also))
        self.wait(12)

        self.play(*[FadeOut(obj) for obj in self.mobjects])


def get_bell_curve_scene(scene, delay):
    ''' Graph of the PDF of the normal distribution '''
    axes = Axes(x_range=(-4, 4, 1), y_range=(0, 0.5, 0.1),
                axis_config={"include_numbers": True})
    scene.play(Create(axes))
    scene.wait(1)

    bell_curve = axes.plot(
        lambda x: math.exp(-x ** 2 / 2) / math.sqrt(2 * PI), x_range=(-4, 4, 0.05))
    bell_curve.set_color(YELLOW)
    scene.play(Create(bell_curve))
    scene.wait(delay)

    scene.play(*[FadeOut(obj) for obj in scene.mobjects])


class BellCurve1(Scene):
    def construct(self):
        get_bell_curve_scene(self, 34)


class RotationInvariant(Scene):
    def construct(self):
        # Proof that the standard multivariate Gaussian is rotationally invariant
        q_mat = MathTex(
            r"Q = \begin{bmatrix} q_{1, 1} & \ldots & q_{1, n} \\ \vdots & \ddots & \vdots \\ q_{n, 1} & \ldots & q_{n, n} \end{bmatrix}")
        x_vec = MathTex(
            r"X = \begin{bmatrix} X_1 \\ \vdots \\ X_n \end{bmatrix}, X_i \sim \mathcal{N}(0, 1)")
        y_qx = MathTex(
            r"Y = \begin{bmatrix} Y_1 \\ \vdots \\ Y_n \end{bmatrix} = QX = \begin{bmatrix} q_{1, 1} & \ldots & q_{1, n} \\ \vdots & \ddots & \vdots \\ q_{n, 1} & \ldots & q_{n, n} \end{bmatrix} \begin{bmatrix} X_1 \\ \vdots \\ X_n \end{bmatrix}, X_i \sim \mathcal{N}(0, 1)")
        y1 = MathTex(r"Y_1 = q_{1, 1} X_1 + \ldots + q_{1, n} X_n")
        q_sq = MathTex(r"q_{1, 1}^2 + \ldots + q_{1, n}^2 = 1")
        mu_y1 = MathTex(
            r"E(Y_1) = q_{1,1} E(X_1) + \ldots + q_{1,n} E(X_n) = 0")
        var_y1 = MathTex(
            r"\mathrm{Var}(Y_1) = q_{1,1}^2 \mathrm{Var}(X_1) + \ldots + q_{1,n}^2 \mathrm{Var}(X_n) = q_{1,1}^2 + \ldots + q_{1,n}^2 = 1")
        y1_normal = MathTex(r"Y_1 \sim \mathcal{N}(0, 1)")
        cov = MathTex(r"& \mathrm{Cov}(Y_1, Y_2) \\ & = E[Y_1 Y_2] - E[Y_1] E[Y_2] \\ & = E[(q_{1,1} X_1 + \ldots + q_{1,n} X_n)(q_{2,1} X_1 + \ldots + q_{2,n} X_n)] \\ & = q_{1,1} q_{2,1} E[X_1^2] + \ldots + q_{1,n} q_{2,n} E[X_n^2] \\ & + \sum_{i,j \leq n, i \neq j} q_{1,i} q_{2,j} E[X_i X_j] \\ & = E[X_1^2] (q_{1,1} q_{2,1} + \ldots + q_{1,n} q_{2,n}) \\ & + \sum_{i,j \leq n, i \neq j} q_{1,i} q_{2,j} E[X_i] E[X_j] \\ & = 0")
        cov.scale(0.7)
        y_normal = Text("Y = QX is still a standard multivariate normal",
                        t2s={"Y = QX": ITALIC})
        y_normal.scale(0.7)  # Text objects are bigger than MathTex objects

        qx_group = VGroup(q_mat, x_vec)
        qx_group.arrange(DOWN)
        self.play(FadeIn(qx_group))
        self.wait(13)

        self.play(ReplacementTransform(q_mat, y_qx), FadeOut(x_vec))
        self.wait(4)

        y1.next_to(y_qx, DOWN)
        self.play(y_qx.animate.shift(UP), FadeIn(y1))
        self.wait(5)

        self.play(FadeOut(y_qx), y1.animate.shift(UP * 3))
        self.wait(3)

        mu_var_group = VGroup(q_sq, mu_y1, var_y1)
        mu_var_group.scale(0.7)
        mu_var_group.arrange(DOWN)
        mu_var_group.next_to(y1, DOWN)
        self.play(FadeIn(q_sq))
        self.wait(14)

        self.play(FadeIn(mu_y1))
        self.wait(5)

        self.play(FadeIn(var_y1))
        self.wait(2)

        self.play(ReplacementTransform(y1, y1_normal), FadeOut(mu_var_group))
        self.wait(9)

        self.play(y1_normal.animate.to_edge(UP))
        self.wait(1)

        self.play(Write(cov))
        self.wait(1)

        self.play(FadeOut(y1_normal), FadeOut(cov), FadeIn(y_normal))
        self.wait(5)

        # Show image of the joint density function of the 2D Multivariate Normal
        joint_df = ImageMobject("images/Multivariate Normal PDF.png")
        joint_df.scale(1.2)
        self.play(FadeOut(y_normal), FadeIn(joint_df))
        self.wait(12)

        self.play(*[FadeOut(obj) for obj in self.mobjects])


class UnifProof(Scene):
    def construct(self):
        r = 3
        circle = Circle(radius=r)
        self.play(Create(circle))
        self.wait(1)

        vec1 = play_normalize_vec(self, [-1.44 * r, 0.67 * r, 0], r)
        vec2 = play_normalize_vec(self, [-0.36 * r, -0.29 * r, 0], r)
        self.wait(6)
        self.play(FadeOut(vec1), FadeOut(vec2))

        eps = 0.1
        diag_arc = get_diag_arc(r, eps, GREEN_D)
        diag_arc_copy = diag_arc.copy()
        a_label = MathTex(r"A")
        a_label.move_to(diag_arc.get_center() + UR * MED_LARGE_BUFF)

        self.play(Create(diag_arc), run_time=0.5)
        self.play(Indicate(diag_arc), run_time=0.5)
        self.play(Write(a_label), run_time=0.5)
        self.wait(0.5)

        vert_arc = get_vert_arc(r, eps, BLUE_D)
        b_label = MathTex(r"B")
        b_label.next_to(vert_arc, UP)

        self.play(Create(vert_arc), run_time=0.5)
        self.play(Indicate(vert_arc), run_time=0.5)
        self.play(Write(b_label), run_time=0.5)
        self.wait(5)

        self.play(FadeOut(vert_arc))
        self.play(Rotate(diag_arc, PI / 4, about_point=ORIGIN))
        self.wait(17)

        self.play(FadeIn(diag_arc_copy))
        self.wait(18)

        self.play(*[FadeOut(obj) for obj in self.mobjects])


# Play BellCurve scene again to mention Box-Muller
class BellCurve2(Scene):
    def construct(self):
        get_bell_curve_scene(self, 24)


class ChooseDist(Scene):
    def construct(self):
        # Title
        title = Text("Cumulative Distribution Function (CDF)")
        title.to_edge(UP)
        self.play(Write(title))
        self.wait(19)

        # F_X(r)
        fx_1 = MathTex(r"F_X(r) = P(X \leq r)")
        fx_2 = MathTex(
            r"= \frac{\text{Vol(Ball with radius $r$)}}{\text{Vol(Ball with radius $1$)}}")
        fx_3 = MathTex(r"= r^n")

        fx_group = VGroup(fx_1, fx_2, fx_3)
        fx_group.arrange(DOWN, aligned_edge=LEFT)
        fx_group.to_edge(LEFT)

        self.play(FadeIn(fx_1))
        self.wait(15)

        self.play(FadeIn(fx_2))
        self.wait(9)

        self.play(FadeIn(fx_3))

        # Graph
        axes = Axes(x_range=(0, 1.2, 0.5), y_range=(0, 1.2, 0.5),
                    axis_config={"include_numbers": True},
                    x_length=config.frame_width / 2 - 1,
                    y_length=config.frame_height - 3)

        cdf = axes.plot(
            lambda x: x ** 3, x_range=(0, 1, 0.02))
        cdf.set_color(YELLOW)

        x = 0.9
        y = x ** 3
        dot = Dot(point=axes.coords_to_point(x, y))

        cdf2 = axes.plot(
            lambda x: x ** 10, x_range=(0, 1, 0.02))
        cdf2.set_color(YELLOW)

        graph_group = VGroup(axes, cdf, cdf2, dot)
        graph_group.next_to(fx_group, RIGHT, buff=LARGE_BUFF)
        graph_group.shift(DOWN)

        self.play(Create(axes))
        self.play(Create(cdf))
        self.wait(5)

        self.play(Create(dot))
        self.wait(1)

        self.play(ReplacementTransform(cdf, cdf2),
                  dot.animate.move_to(axes.coords_to_point(x, x ** 10)))
        self.wait(12)

        # Example for volume concentration near surface
        concen_surf = MathTex(
            r"\text{When } n = 500, \\ P(X \geq 0.99) & = 1 - 0.99^{500} \\ & \approx 99.3\%")
        concen_surf.next_to(fx_group, DOWN)
        concen_surf.scale(0.7)
        self.play(FadeIn(concen_surf))
        self.wait(12)

        self.play(FadeOut(concen_surf))
        self.wait(6)

        r_expr = MathTex(r"X = F_X^{-1}(U) = \sqrt[n]{U}")
        r_expr.next_to(fx_group, DOWN)
        self.play(FadeIn(r_expr))
        self.wait(7)

        self.play(*[FadeOut(obj) for obj in self.mobjects])


def random_unit_vector(n):
    ''' Returns a unit vector in R^n chosen uniformly at random. '''
    vec = np.empty(n)
    for i in range(n):
        vec[i] = np.random.normal(0, 1)
    length = np.linalg.norm(vec)
    return vec / length


def normal_sampling(n):
    ''' Returns a random point in an n-dimensional ball. Runs in polynomial time in n. '''
    vec = random_unit_vector(n)
    r = (random.random()) ** (1.0 / n)
    return vec * r


class NormalSampSim2D(Scene):
    def construct(self):

        ns_code_raw = '''
        def random_unit_vector(n):
            vec = np.empty(n)
            for i in range(n):
                vec[i] = np.random.normal(0, 1)
            length = np.linalg.norm(vec)
            return vec / length
        
        def normal_sampling(n):
            vec = random_unit_vector(n)
            r = (random.random()) ** (1.0 / n)
            return vec * r
        '''
        ns_code_rendered = Code(
            code=ns_code_raw,
            language="Python",
            font="Monospace",
            background_stroke_width=0,
            insert_line_no=False)
        self.play(FadeIn(ns_code_rendered))
        self.wait(40)
        self.play(FadeOut(ns_code_rendered))
        self.wait(1)

        r = 3
        p = 3141  # Number of points
        points = np.empty((p, 3))
        for i in range(p):
            points[i] = np.append(r * normal_sampling(2), 0)

        circle = Circle(radius=r)
        self.play(Create(circle))

        dots = VGroup()

        for point in points:
            dots.add(Dot(point, color=BLUE,
                     radius=DEFAULT_DOT_RADIUS / 3))

        self.play(Create(dots), lag_ratio=0.1, run_time=5)
        self.wait(1)

        self.play(*[FadeOut(obj) for obj in self.mobjects])


class NormalSampSim3D(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=PI / 3, theta=PI / 4)
        self.begin_ambient_camera_rotation()

        # axes = ThreeDAxes()
        # self.add(axes)

        r = 3
        p = 3141  # Number of points
        points = np.empty((p, 3))
        for i in range(p):
            points[i] = r * normal_sampling(3)

        sphere = Sphere(radius=r)
        sphere.set_color(GRAY)
        sphere.set_opacity(0.2)
        self.play(Create(sphere))

        dots = VGroup()

        for point in points:
            dots.add(Dot(point, color=BLUE,
                     radius=DEFAULT_DOT_RADIUS / 3))

        self.play(Create(dots), lag_ratio=0.1, run_time=5)
        self.wait(1)

        self.play(*[FadeOut(obj) for obj in self.mobjects])


def fade_and_wait_blist(scene, blist, delays):
    ''' Fade in a bulleted list point by point,
        where the list delays specifies the time to wait after each point.
        Requires: blist and delays have the same length. '''
    for i in range(len(delays)):
        scene.play(FadeIn(blist[i]))
        scene.wait(delays[i])


class Summary(Scene):
    def construct(self):
        title = Text("Recap")
        title.to_edge(UP)
        self.play(Write(title))

        blist = BulletedList("Rejection sampling - huge number of trials",
                             "Spherical coordinates - complicated distribution",
                             "Random direction + random distance",
                             "Random point on sphere - standard multivariate normal",
                             "Random distance from origin - inverse transform sampling")

        delays = [19, 7, 5, 21, 10]

        fade_and_wait_blist(self, blist, delays)

        self.play(*[FadeOut(obj) for obj in self.mobjects])


class Reflection(Scene):
    def construct(self):
        title = Text("Reflection")
        title.to_edge(UP)
        self.play(Write(title))

        blist = BulletedList("Higher-dimensional geometry is interesting",
                             "Making this video is a rewarding experience",
                             "Teaching myself Manim - neat style, cool animations")

        delays = [16, 12, 26]

        fade_and_wait_blist(self, blist, delays)

        self.play(*[FadeOut(obj) for obj in self.mobjects])


class References(Scene):
    def construct(self):
        ref_txt = Text("References")
        ref_txt.to_edge(UP)
        self.play(Write(ref_txt))

        refs = BulletedList(r"Justin's video ``The BEST Way to Find a Random Point in a Circle''",
                            r"\textit{Vector Calculus, Linear Algebra, and Differential Forms: A Unified Approach} by John Hubbard and Barbara Burke Hubbard",
                            r"\textit{Foundations of Data Science} by Avrim Blum, John Hopcroft, and Ravindran Kannan",
                            r"\textit{Introduction to Probability} by David F. Anderson, Timo Sepp\"{a}l\"{a}inen, and Benedek Valk\'{o}",
                            r"This website helped me generate the graph of the PDF of the standard multivariate normal: https://scipython.com/blog/visualizing-the-bivariate-gaussian-distribution/",
                            r"Manim, numpy, and matplotlib documentations and tutorials")
        refs.scale_to_fit_width(12)
        refs.next_to(ref_txt, DOWN, buff=MED_LARGE_BUFF)
        self.play(Write(refs, run_time=5))
        self.wait(3)
        self.play(*[FadeOut(obj) for obj in self.mobjects])
