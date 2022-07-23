from telnetlib import DO
from turtle import down
from manim import *


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
