import random
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


class BallVolume(Scene):
    def construct(self):
        # self.add(NumberPlane())  # Debug

        t_int = MathTex(r"c_n = \int_{-1}^1 (1 - t^2)^{\frac{n-1}{2}} \, dt")
        self.play(Create(t_int))
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
