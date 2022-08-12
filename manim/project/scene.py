import math
import random
import numpy as np
from manim import *


class Intro(Scene):
    def construct(self):
        hi = Text("Hey everyone, I'm James.")
        justin = ImageMobject("images/Justin's Video.png")
        hi.to_edge(UP, buff=0.2)
        justin.scale(0.9)
        justin.next_to(hi, DOWN)
        self.play(Create(hi))
        self.wait(9)

        self.play(FadeIn(justin))
        self.wait(15)

        self.play(FadeOut(hi), FadeOut(justin))
        self.wait(1)

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

        self.play(Create(problem))
        self.wait(15)

        self.play(Create(circ_def))
        self.wait(7)

        self.play(Create(ball_def))
        self.wait(10)

        # Fade out everything
        self.play(*[FadeOut(obj) for obj in self.mobjects])


class RejSampling(Scene):
    def construct(self):
        rej_samp_title = Text("Rejection Sampling")
        rej_samp_title.to_edge(UP, buff=0.2)
        self.play(Create(rej_samp_title))
        self.wait(2)

        r = 3
        square = Square(side_length=r * 2)
        circle = Circle(radius=r)

        self.play(Create(square), Create(circle))
        self.wait(1)

        num_points = 5  # Must be greater than 1
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
        # Fade out all objects
        self.play(*[FadeOut(obj) for obj in self.mobjects])
        self.wait(1)

        rs_code_raw = '''
        def rejection_sampling(n):
            trials = 0
            while True:
                trials += 1
                vec = np.empty(n)
                for i in range(n):
                    vec[i] = random.random()
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
        self.wait(2)

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
        self.wait(2)

        self.play(FadeIn(trials_log))

        # Fade out everything
        self.play(*[FadeOut(obj) for obj in self.mobjects])


class BallVolume(Scene):
    def construct(self):
        # self.add(NumberPlane())  # Debug

        r_x1 = MathTex(r"r = \sqrt{1 - x_1^2}")
        beta_n = MathTex(r"\beta_n &= \int_{-1}^1 r^{n-1} \beta_{n-1} \, dx_1")
        beta_n_2 = MathTex(
            r"\beta_n &= \beta_{n-1} \int_{-1}^1 (1 - x_1^2)^{\frac{n-1}{2}} \, dx_1")
        r_x1.next_to(beta_n, UP)

        self.play(FadeIn(r_x1))
        self.play(Create(beta_n))
        self.wait(1)

        self.play(ReplacementTransform(beta_n, beta_n_2), FadeOut(r_x1))
        self.wait(1)

        t_int = MathTex(r"c_n = \int_{-1}^1 (1 - t^2)^{\frac{n-1}{2}} \, dt")
        t_int.save_state()
        t_int.next_to(beta_n_2, DOWN)
        t_int.scale(0.7)

        beta_cn = MathTex(r"\beta_n = c_n \beta_{n - 1}")
        beta_cn.save_state()
        beta_cn.next_to(t_int, DOWN)
        beta_cn.scale(0.7)

        self.play(FadeIn(t_int))
        self.play(FadeIn(beta_cn))
        self.wait(1)

        self.play(FadeOut(beta_n_2),
                  FadeOut(beta_cn),
                  Restore(t_int))
        self.wait(1)

        trig_sub = MathTex(r"t = \sin(\theta)")
        trig_sub.next_to(t_int, DOWN)
        trig_sub.scale(0.7)
        dt = MathTex(r"dt = \cos(\theta) \, d\theta")
        dt.next_to(trig_sub, DOWN)
        dt.scale(0.7)
        self.play(FadeIn(trig_sub))
        self.play(FadeIn(dt))
        self.wait(1)

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
        self.wait(1)
        self.play(FadeOut(trig_int_n_2))
        self.play(trig_int_n.animate.shift(DOWN))

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
        self.wait(1)

        ibp = MathTex(
            r"c_n = \cos^{n-1}(\theta) \sin(\theta) \bigg|_{-\frac{\pi}{2}}^{\frac{\pi}{2}} + \int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} \sin^2(\theta) (n-1) \cos^{n-2}(\theta) \, d\theta")
        ibp_text = Text("Integration by Parts", font_size=36)
        ibp_text.next_to(ibp, UP)
        self.play(FadeOut(u_dv), FadeOut(du_v),
                  ReplacementTransform(trig_int_n, ibp), FadeIn(ibp_text))
        self.wait(1)

        cos_is_zero = MathTex(
            r"\cos \left(-\frac{\pi}{2} \right) = \cos{\frac{\pi}{2}} = 0")
        cos_is_zero.next_to(ibp, DOWN)
        cos_is_zero.scale(0.7)
        self.play(FadeOut(ibp_text))
        self.play(FadeIn(cos_is_zero))
        self.wait(1)

        ibp_one_term = MathTex(
            r"c_n = (n-1) \int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} \sin^2(\theta) \cos^{n-2}(\theta) \, d\theta")
        self.play(ReplacementTransform(
            ibp, ibp_one_term), FadeOut(cos_is_zero))
        self.wait(1)

        sin_cos_id = MathTex(r"\sin^2(\theta) = 1 - \cos^2(\theta)")
        sin_cos_id.next_to(ibp_one_term, DOWN)
        sin_cos_id.scale(0.7)
        self.play(FadeIn(sin_cos_id))
        self.wait(1)

        ibp_split = MathTex(
            r"c_n = (n-1) \left[\int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} \cos^{n-2}(\theta) \, d\theta - \int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} \cos^n(\theta) \, d\theta\right]")
        self.play(ReplacementTransform(
            ibp_one_term, ibp_split), FadeOut(sin_cos_id))
        self.wait(1)

        self.play(ibp_split.animate.shift(UP))
        trig_int_n = MathTex(
            r"c_n &= \int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} \cos^{n}(\theta) \, d\theta")
        trig_int_n.next_to(ibp_split, DOWN)
        trig_int_n_2.next_to(trig_int_n, DOWN)
        self.play(FadeIn(trig_int_n))
        self.play(FadeIn(trig_int_n_2))
        self.wait(1)

        ibp_rearrange = MathTex(r"c_n = (n-1) c_{n-2} - (n-1) c_n")
        self.play(ReplacementTransform(
            ibp_split, ibp_rearrange), FadeOut(trig_int_n), FadeOut(trig_int_n_2))
        self.wait(1)

        result = MathTex(r"c_n = \frac{n-1}{n}c_{n-2}")
        self.play(ReplacementTransform(ibp_rearrange, result))
        self.wait(1)

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
        self.wait(1)

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
        self.wait(1)

        kappa_rel = MathTex(
            r"\kappa_{n - 1} \stackrel{\times 2}{\longrightarrow} \kappa_n")
        kappa_rel.next_to(beta_cn, DOWN)
        self.play(FadeIn(kappa_rel))
        self.wait(1)

        beta_rel = MathTex(
            r"\beta_{n - 1} \stackrel{\times c_n}{\longrightarrow} \beta_n")
        beta_rel.next_to(kappa_rel, DOWN, buff=MED_LARGE_BUFF)
        self.play(FadeIn(beta_rel))

        bk_group = VGroup(kappa_rel, beta_rel)

        ratio_rel = MathTex(
            r"\frac{\beta_{n - 1}}{\kappa_{n-1}} \stackrel{\times \frac{c_n}{2}}{\longrightarrow} \frac{\beta_n}{\kappa_n}"
        )
        ratio_rel.next_to(beta_cn, DOWN, buff=MED_LARGE_BUFF)
        self.play(FadeOut(bk_group),
                  FadeIn(ratio_rel))
        self.wait(1)

        prob = MathTex(
            r"P(\text{Random point in ball}) = \frac{\mathrm{Vol(Ball)}}{\mathrm{Vol(Box)}} = \frac{\beta_n}{\kappa_n}")
        prob.next_to(ratio_rel, DOWN)
        self.play(FadeIn(prob))
        self.wait(1)

        ex = MathTex(r"E[X] = \frac{1}{P(\text{Random point in ball)}}")
        incr_fast = Text("increases faster than exponentially", color=RED)
        ex_group = VGroup(ex, incr_fast)
        ex_group.arrange(DOWN)
        ex_group.next_to(prob, DOWN)
        self.play(FadeIn(ex))
        self.play(GrowFromCenter(incr_fast))

        # Fade out everything
        self.play(*[FadeOut(obj) for obj in self.mobjects])


class SphereDistIntro(Scene):
    def construct(self):
        r = 3
        circle = Circle(radius=r)
        self.play(Create(circle))

        theta_0 = PI / 3
        unit_vec = Vector([r * math.cos(theta_0), r * math.sin(theta_0), 0])
        self.play(Create(unit_vec))
        self.play(Rotate(unit_vec,
                         angle=2 * PI,
                         about_point=ORIGIN,
                         run_time=2,
                         rate_func=smooth))
        self.wait(1)
        self.play(ApplyMethod(unit_vec.put_start_and_end_on,
                              [0, 0, 0],
                              [0.2 * r * math.cos(theta_0),
                               0.2 * r * math.sin(theta_0), 0],
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
        self.wait(1)

        self.play(FadeIn(sphere_def))
        self.wait(1)

        # Fade out everything
        self.play(*[FadeOut(obj) for obj in self.mobjects])


class PointOnSphereWrong(Scene):
    def construct(self):
        r = 3
        square = Square(side_length=r * 2)
        circle = Circle(radius=r)

        self.play(Create(square), Create(circle))
        self.wait(1)

        def play_normalize_vec(arr):
            ''' Creates a vector based on arr and then normalizes it.
                Returns the unit vector created.
            '''
            unit_arr = arr / np.linalg.norm(arr) * r
            vec = Vector(arr)
            unit_vec = Vector(unit_arr)
            self.play(Create(vec))
            self.wait(1)
            self.play(ReplacementTransform(vec, unit_vec))
            self.wait(1)
            return unit_vec  # So that we can fade it out later

        arrs = [np.array([0.6 * r, -0.4 * r]),
                np.array([-0.93 * r, -0.78 * r])]

        unit_vecs = []

        for arr in arrs:
            unit_vecs.append(play_normalize_vec(arr))

        self.play(*[FadeOut(unit_vec) for unit_vec in unit_vecs])

        diag = Line(ORIGIN, [r, r, 0])
        vert = Line(ORIGIN, [0, r, 0])
        diag_label = MathTex(r"\sqrt{2}")
        vert_label = MathTex(r"1")
        diag_label.move_to(diag.get_center() + DR * MED_SMALL_BUFF)
        vert_label.next_to(vert, buff=SMALL_BUFF)
        diag_group = VGroup(diag, diag_label)
        vert_group = VGroup(vert, vert_label)

        self.play(Create(diag_group))
        self.wait(1)
        self.play(Create(vert_group))
        self.wait(1)
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

        diag_arc_start = math.atan2(1 - 2 * eps, 1)
        diag_arc_end = math.atan2(1, 1 - 2 * eps)
        diag_arc = Arc(arc_center=ORIGIN,
                       radius=r,
                       start_angle=diag_arc_start,
                       angle=diag_arc_end - diag_arc_start)
        diag_arc.set_color(GREEN_D)

        vert_arc_start = math.atan2(1, eps)
        vert_arc_end = math.atan2(1, -eps)
        vert_arc = Arc(arc_center=ORIGIN,
                       radius=r,
                       start_angle=vert_arc_start,
                       angle=vert_arc_end - vert_arc_start)
        vert_arc.set_color(BLUE_D)

        diag_arr = np.array([0.87 * r, (0.87 - 0.2 * eps) * r, 0])
        play_normalize_vec(diag_arr)

        self.play(Create(diag_arc))
        self.play(Indicate(diag_arc))
        self.wait(1)

        vert_arr = np.array([-eps * 0.05 * r, 0.47 * r, 0])
        play_normalize_vec(vert_arr)

        self.play(Create(vert_arc))
        self.play(Indicate(vert_arc))
        self.wait(1)

        graph_group = Group(*[obj for obj in self.mobjects])
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
        txt_group.arrange(DOWN)
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
        self.wait(1)

        self.play(FadeIn(root_2_sqr))
        self.wait(1)

        self.play(ReplacementTransform(root_2_sqr, root_n_pow),
                  ReplacementTransform(prob_eq_area, prob_eq_vol))
        self.wait(1)

        self.play(GrowFromCenter(incr_fast))
        self.wait(1)

        # Fade out everything
        self.play(*[FadeOut(obj) for obj in self.mobjects])


class RotateSquare(Scene):
    def construct(self):
        r = 2
        square = Square(side_length=r * 2)
        self.play(Create(square))
        self.play(Rotate(square, angle=PI / 4))
        self.wait(1)
        self.play(Uncreate(square))


class RotatePlane(LinearTransformationScene):
    def __init__(self):
        LinearTransformationScene.__init__(
            self
        )

    def construct(self):
        q_mat = MathTex(
            r"Q = \begin{bmatrix} \cos(\theta) & -\sin(\theta) \\ \sin(\theta) & \cos(\theta) \end{bmatrix}")

        q_mat.to_corner(UR)
        self.add_foreground_mobject(q_mat)  # Foreground mobjects don't move

        theta = PI / 3
        matrix = [[math.cos(theta), -math.sin(theta)],
                  [math.sin(theta), math.cos(theta)]]
        self.apply_matrix(matrix)
        self.wait()


class PointOnSphereRight(Scene):
    def construct(self):
        q_mat = MathTex(
            r"Q = \begin{bmatrix} \cos(\theta) & -\sin(\theta) \\ \sin(\theta) & \cos(\theta) \end{bmatrix}")

        orthonormal = MathTex(
            r"\begin{bmatrix} \cos(\theta) \\ \sin(\theta) \end{bmatrix} \text{ and } \begin{bmatrix} -\sin(\theta) \\ \cos(\theta) \end{bmatrix} \text{ are orthonormal}")
        norm_1 = MathTex(r"\sqrt{\cos^2(\theta) + \sin^2(\theta)} = 1")
        dot_0 = MathTex(
            r"\begin{bmatrix} \cos(\theta) \\ \sin(\theta) \end{bmatrix} \cdot \begin{bmatrix} -\sin(\theta) \\ \cos(\theta) \end{bmatrix} = 0")

        rows_also = Text(
            "Similarly, the row vectors of Q are also orthonormal.", t2s={"Q": ITALIC})
        rows_also.scale(0.6)

        q_group = VGroup(q_mat, orthonormal, norm_1, dot_0, rows_also)
        q_group.arrange(DOWN)

        self.play(FadeIn(q_mat))
        self.wait(1)

        self.play(FadeIn(orthonormal), FadeIn(norm_1), FadeIn(dot_0))
        self.wait(1)

        self.play(FadeIn(rows_also))
        self.wait(1)


def random_unit_vector(n):
    ''' Returns a unit vector in R^n chosen uniformly at random. '''
    vec = np.empty(n)
    for i in range(n):
        vec[i] = np.random.normal(0., 1.)
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
                vec[i] = np.random.normal(0., 1.)
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
        self.wait(2)
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

        self.play(Create(dots), lag_ratio=0.1, run_time=2)
        self.wait(1)


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

        self.play(Create(dots), lag_ratio=0.1, run_time=10)
        self.wait(1)


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
                            r"Manim, numpy, and matplotlib documentations")
        refs.scale_to_fit_width(12)
        refs.next_to(ref_txt, DOWN)
        self.play(Write(refs, run_time=5))
        self.wait(5)
