import math
from manim import *
import numpy as np

# As shown in https://www.youtube.com/watch?v=SYGLhIKoDc8

# Quality settings - easily adjust for draft vs final rendering
class QualitySettings:
    # Set DRAFT to True for faster rendering during development
    # Set to False for final high-quality output
    DRAFT = False
    
    @classmethod
    def get_settings(cls):
        if cls.DRAFT:
            return {
                "cylinder_resolution": (10, 10),      # Low resolution for faster rendering
                "surface_resolution": (20, 20),     # Lower resolution surface
                "curve_sample_points": 50,          # Fewer points for curves
                "refinement_steps": 1,              # Fewer refinement iterations
                "fill_opacity": 0.6,                # Less opacity for better visibility in draft
                "run_time_scale": 0.5,              # Faster animations for quicker previewing
            }
        else:
            return {
                "cylinder_resolution": (20, 40),  # High resolution for final output
                "surface_resolution": (30, 30),   # Detailed surface
                "curve_sample_points": 200,         # Smooth curves
                "refinement_steps": 2,              # More refinements for better transition
                "fill_opacity": 0.7,                # Full opacity for final render
                "run_time_scale": 1.5,              # Normal animation speed
            }

class VolumeCalculationScene(ThreeDScene):
    def construct(self):
        # Get quality settings at the beginning
        settings = QualitySettings.get_settings()
        
        thin_height = 0.3  # Very thin initial height
        dx = thin_height
        num_slices = int(2/dx)
        x_values = np.linspace(-0.8, 0.8, num_slices)
        x_circle = x_values[3]  # Choose a specific x value for the circle
        # Define the vase profile function once at the top
        def vase_function(x):
            return 0.2+0.1*(x+1)+0.8*math.sin(0.6*PI*(x+1))**2
        
        # Start with camera already in position to see the z-y plane
        self.set_camera_orientation(phi=-90 * DEGREES, theta=0, gamma=90 * DEGREES)
        
        # --- Szene 1: Kreis in der z-y-Ebene zeichnen ---
        # Create circle based on axis scaling to ensure proper size
        # First create the 3D axes early but don't show them yet
        axes_3d = ThreeDAxes(
            x_range=[-1.5, 1.5, 0.5],
            y_range=[-1.5, 1.5, 0.5],
            z_range=[-1.5, 1.5, 0.5],
            axis_config={"include_tip": True, "include_ticks": False},
            x_length=5, y_length=5, z_length=5
        )
        axes_3d.shift(DOWN)  # Shift axes down for more space at the top

        # Add properly oriented labels AFTER shifting the axes
        x_label = MathTex("x").scale(0.7).next_to(axes_3d.get_x_axis(), RIGHT)
        y_label = MathTex("y").scale(0.7).next_to(axes_3d.get_y_axis(), UP)
        z_label = MathTex("z").scale(0.7).next_to(axes_3d.get_z_axis(), OUT)

        # Rotate labels to face camera
        x_label.rotate(-90 * DEGREES, axis=UP)
        y_label.rotate(-90 * DEGREES, axis=UP)
        z_label.rotate(-90 * DEGREES, axis=UP)

        radius = vase_function(x_circle)  # Radius based on our function

        # Create circle in 3D coordinates with quality-adjusted resolution
        circle = ParametricFunction(
            lambda t: axes_3d.c2p(x_circle, radius*np.cos(t), radius*np.sin(t)),
            t_range=(0, TAU, TAU/settings["curve_sample_points"]),  # Adjust point density
            color=WHITE
        )

        # Create dot on circle
        dot = Dot(color=YELLOW)
        dot.move_to(circle.point_from_proportion(0))
        dot.rotate(-90 * DEGREES, axis=UP)  # Rotate to face camera

        # Create a STATIC trace first - this will persist
        static_trace = circle.copy().set_stroke(color=RED, width=2, opacity=0.8)

        # Create the dynamic trace for animation
        trace = TracedPath(dot.get_center, stroke_color=RED, stroke_width=2)
        self.add(trace, dot)

        # Create updater for dot to follow circle
        def update_dot(d, alpha):
            d.move_to(circle.point_from_proportion(alpha))
            return d

        # Animate dot with quality-adjusted runtime
        self.play(
            UpdateFromAlphaFunc(dot, update_dot), 
            run_time=3 * settings["run_time_scale"]
        )
        self.wait(0.5 * settings["run_time_scale"])

        # Once the animation is done, replace the dynamic trace with the static one
        self.remove(trace)
        self.add(static_trace)

        # Zeige den Text "Fläche = π r²" - with consistent height
        area_text = MathTex(r"\text{Fläche} = \pi r^2")
        area_text.to_edge(UP).shift(DOWN * 0.5)  # Add the DOWN shift for consistency
        area_text.rotate(-90 * DEGREES, axis=UP)  # Rotate to face camera
        # Show area formula with quality-adjusted runtime
        self.play(Write(area_text), run_time=1 * settings["run_time_scale"])
        self.wait(1 * settings["run_time_scale"])

        # --- Scene 1: Create a very thin cylinder (representing a circular slice) ---

        # Create a cylinder at position z_circle
        # The cylinder extends along the x-axis and its circular face is in the z-y plane
        cylinder_slice = Cylinder(
            radius=radius*1.66666666666667,  # Adjust for aspect ratio
            height=thin_height,
            direction=RIGHT,  # Extends along x-axis
            fill_opacity=settings["fill_opacity"],
            fill_color=RED,
            stroke_width=0,
            checkerboard_colors=[RED, RED_D],
            resolution=settings["cylinder_resolution"]
        )
        cylinder_slice.move_to(axes_3d.c2p(x_circle, 0, 0))

        # Show the coordinate system with quality-adjusted runtime
        axes_labels = VGroup(x_label, y_label)
        self.play(
            Create(axes_3d), 
            Write(axes_labels), 
            run_time=1 * settings["run_time_scale"]
        )
        self.wait(1 * settings["run_time_scale"])

        # Prepare the new text but don't display it yet - with consistent height
        volume_text_init = MathTex(r"\text{Volumen} = \pi r^2\, \Delta x")
        volume_text_init.to_edge(UP).shift(DOWN * 0.5)  # Add the DOWN shift for consistency

        def rotate_system(degrees):
            # First move the camera
            anims = [
                Rotate(x_label, angle=degrees * DEGREES, axis=UP),
                Rotate(y_label, angle=degrees * DEGREES, axis=UP),
                Rotate(dot, angle=degrees * DEGREES, axis=UP),
                Transform(static_trace, cylinder_slice),
                Transform(area_text, volume_text_init)  # Transform old text to new text
            ]
            # Camera rotation with quality-adjusted runtime
            self.move_camera(
                phi=(-90+degrees)*DEGREES, 
                theta=0, 
                gamma=degrees * DEGREES, 
                run_time=2 * settings["run_time_scale"], 
                added_anims=anims
            )
            self.wait(1 * settings["run_time_scale"])

        rotate_system(90)
        self.wait(1)
        # --- Szene 2: Funktion f(x) einblenden ---
        # Create 2D axes that precisely match the 3D axes scale after rotation
        # axes = Axes(
        #     x_range=[-1.5, 1.5, 0.5],  # Match the 3D axes range exactly
        #     y_range=[-1.5, 1.5, 0.5],  # Match the 3D axes range exactly
        #     x_length=5,                # Match the 3D axes length exactly
        #     y_length=5,                # Match the 3D axes length exactly
        #     tips=True,
        #     axis_config={"include_numbers": True}
        # )

        # Plot the vase function using our defined function
        func = axes_3d.plot(
            lambda x: vase_function(x), 
            color=BLUE,
            x_range=[-1, 1, 2/settings["curve_sample_points"]]  # Adjust curve resolution
        )
        minus_func = axes_3d.plot(
            lambda x: -vase_function(x), 
            color=BLUE,
            x_range=[-1, 1, 2/settings["curve_sample_points"]]  # Adjust curve resolution
        )
        func_label = MathTex("f(x)")  # Fix label to match variable
        func_label.next_to(func, RIGHT)

        # Actually display the axes, function and label
        # self.play(FadeIn(axes))
        # Show function with quality-adjusted runtime
        self.play(
            Create(func), 
            Create(minus_func), 
            Write(func_label), 
            run_time=1 * settings["run_time_scale"]
        )
        self.wait(1 * settings["run_time_scale"])

        # --- Szene 3: Volumenscheibe ---
        dx = thin_height
        self.play(FadeOut(dot))

        # Update text - transform from the previous state - with consistent height
        volume_text = MathTex(r"\text{Volumen} = \pi f(x)^2\, \Delta x")
        volume_text.to_edge(UP).shift(DOWN * 0.5)  # Add the DOWN shift for consistency
        # Update text with quality-adjusted runtime
        self.play(
            Transform(area_text, volume_text), 
            run_time=1 * settings["run_time_scale"]
        )
        self.wait(1 * settings["run_time_scale"])

        # --- Szene 4: Summe von Volumenscheiben ---
        # Create multiple cylinder slices along the x-axis
        multi_slices = VGroup()

        # Define a consistent color scheme for ALL cylinders throughout the animation
        colors = [YELLOW, GREEN, PURPLE,RED, BLUE, ORANGE]

        # Find the index closest to our original cylinder position
        closest_idx = np.argmin(np.abs(x_values - x_circle))

        # Create completely new cylinders for ALL positions (including the original one)
        for i, x in enumerate(x_values):
            # Create a new cylinder for every position - no special handling for original
            radius = vase_function(x) * 1.66666666666667
            c = colors[i % len(colors)]  # Use the consistent color scheme
            slice_cylinder = Cylinder(
                radius=radius,
                height=dx,
                direction=RIGHT,
                fill_opacity=settings["fill_opacity"],
                fill_color=c,
                checkerboard_colors=[c, c],
                stroke_width=0,
                resolution=settings["cylinder_resolution"]
            )
            slice_cylinder.move_to(axes_3d.c2p(x, 0, 0))
            multi_slices.add(slice_cylinder)

        # Show the complete set of new cylinders
        self.play(Create(multi_slices), run_time=1 * settings["run_time_scale"])
        # IMPORTANT: Remove the original cylinder_slice completely from the scene
        self.remove(cylinder_slice)
        self.remove(static_trace)  # Also remove the static_trace to be safe

        # First calculate the actual number of slices we're using
        num_slices = len(x_values)

        # SPLIT the formula into two parts: static label and dynamic summation
        # This way only the summation part changes, not the "Volumen ="
        volume_label = MathTex(r"\text{Volumen} =")
        volume_label.to_edge(UP).to_edge(LEFT).shift(RIGHT * 1.5 + DOWN * 0.5)  # Already has DOWN shift

        # Create the summation part separately
        sum_formula = MathTex(r"\sum_{i=1}^{" + str(num_slices) + r"} \pi f(x_i)^2\, \Delta x")
        sum_formula.next_to(volume_label, RIGHT)

        # Update sum text with quality-adjusted runtime - add both parts separately
        self.play(
            Transform(area_text, volume_label),  # Transform the old text to the new label
            FadeIn(sum_formula),  # Fade in the summation part
            run_time=1 * settings["run_time_scale"]
        )
        self.wait(2 * settings["run_time_scale"])

        # --- Final Scene: Limit as dx approaches 0 ---
        # Skip the transition to limit_text with Δx → 0, go directly to refinements
        current_cylinders = multi_slices
        refinement_steps = settings["refinement_steps"]

        for step in range(1, refinement_steps + 1):
            # Double the number of slices each time
            dx_refined = dx / (2**step)
            num_slices_refined = int(2/dx_refined)
            x_values_refined = np.linspace(-0.8, 0.8, num_slices_refined)
            
            # Update ONLY the summation part with the new number of slices
            new_sum_formula = MathTex(
                r"\sum_{i=1}^{" + str(num_slices_refined) + r"} \pi f(x_i)^2\, \Delta x"
            )
            new_sum_formula.next_to(volume_label, RIGHT)  # Position it next to the stable label
            
            refined_cylinders = VGroup()
            
            # Create refinement cylinders with IDENTICAL positioning and coloring logic
            for i, x in enumerate(x_values_refined):
                radius = vase_function(x) * 1.66666666666667
                c = colors[i % len(colors)]
                slice_cylinder = Cylinder(
                    radius=radius,
                    height=dx_refined,
                    direction=RIGHT,
                    fill_opacity=settings["fill_opacity"],
                    fill_color=c,
                    checkerboard_colors=[c, c],
                    stroke_width=0,
                    resolution=settings["cylinder_resolution"]
                )
                slice_cylinder.move_to(axes_3d.c2p(x, 0, 0))
                refined_cylinders.add(slice_cylinder)
            
            # Smoothly transform both the cylinders AND the text
            self.play(
                ReplacementTransform(current_cylinders, refined_cylinders),
                Transform(sum_formula, new_sum_formula),  # Only transform the summation part
                run_time=1.5 * settings["run_time_scale"]
            )
            
            self.wait(1 * settings["run_time_scale"])
            current_cylinders = refined_cylinders

        # When all refinements are done, go to the final form with limit and integral
        # Again split into parts to maintain stable elements
        final_limit_formula = MathTex(
            r"\lim_{n \to \infty} \sum_{i=1}^{n} \pi f(x_i)^2\, \Delta x = \int_{a}^{b} \pi f(x)^2\, dx"
        )
        final_limit_formula.next_to(volume_label, RIGHT)

        # Transform just the sum part to the final form
        self.play(
            Transform(sum_formula, final_limit_formula),
            run_time=1 * settings["run_time_scale"]
        )
        self.wait(1 * settings["run_time_scale"])

        # Create the revolution surface directly using the original axes
        revolution_surface = Surface(
            lambda u, v: axes_3d.c2p(
                u,
                vase_function(u) * np.cos(v),
                vase_function(u) * np.sin(v)
            ),
            u_range=[-1, 1],
            v_range=[0, TAU],
            resolution=settings["surface_resolution"],
            fill_opacity=settings["fill_opacity"],
            fill_color=BLUE_D,
            stroke_width=0,
            stroke_color=BLUE_E,
            checkerboard_colors=[BLUE_D, BLUE_E]
        )

        # Fade out all cylinders
        self.play(FadeOut(current_cylinders), Create(revolution_surface),
            run_time=2 * settings["run_time_scale"])

        # Instead of using the non-existent cache clearing method, do this:
        self.remove(current_cylinders)  # Explicitly remove from scene

        # Also ensure any other potential cylinder references are gone
        self.remove(cylinder_slice)  # Remove original cylinder if it somehow persists
        self.remove(multi_slices)    # Remove original multi_slices

        self.wait(2 * settings["run_time_scale"])

        # Continue with camera rotation
        # Final camera movement to show the 3D shape from different angles
        # Camera rotation with quality-adjusted parameters
        self.begin_ambient_camera_rotation(
            rate=0.2 * settings["run_time_scale"], 
            about="phi"
        )
        self.wait(5 * settings["run_time_scale"])
        self.stop_ambient_camera_rotation()

        # Final fade out
        # Final fadeout with quality-adjusted runtime
        self.play(
            FadeOut(revolution_surface),
            FadeOut(volume_label),      # Fade out the stable label
            FadeOut(sum_formula),       # Fade out the formula part
            FadeOut(area_text),         # Also fade out the area_text which was transformed to volume_label
            FadeOut(axes_3d),
            FadeOut(axes_labels),
            FadeOut(func),
            FadeOut(minus_func),
            FadeOut(func_label),
            run_time=1 * settings["run_time_scale"]
        )
        self.wait(1 * settings["run_time_scale"])