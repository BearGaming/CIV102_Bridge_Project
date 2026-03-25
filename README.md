This project includes a python script that was written for the purpose of rapidly iterating through different kinds of trapezoidal cross sections given constraints on maximum height, width and amount of provided materials. Our team was attracted to the trapezoid because it is a well tested cross section featured in many highway and railway bridges and uses the disproportionate allocation of material near the "top" of the cross section to bring the centroid further up. This is desirable as it better distributed the stress due to the bending moment between the material crushing at the deck near the deck (5MPa) and failing in tension at the bottom (30MPa).

## Example Output
Generated images of generated cross sections can be found [here](./examples)
<img width="640" height="480" alt="Height = 40" src="https://github.com/user-attachments/assets/6ae26216-619c-42a6-a503-34b294dd8590" />
This image shows how the cross section is represented by the program as a collection of parallelograms each with their own centroids, as well as the overal centroid with a dashed line and reported value in the legend (ybar), calculated using the first moment of area and total area.
