# Implicit Surface Model Generator

Generates non-manifold meshes for the triply periodic minimal surfaces Gyroid and Fisher-Koch S.
Developed to enable 3D printing of scaffolds for bone tissue reconstruction.
The project was completed at Colorado State University under the advisement of Professors David Prawel and Clayton Shonkwiler for partial fulfillment of an undergraduate research course.

# Citation

If you find this software insightful and use it in your research or elsewhere, please cite it.
Github has a new feature to cite a repository that you should see in the top right of the webpage.

```
@software{Gunn_Implicit_Surface_Model_2021,
  author = {Gunn, Erin},
  month = {12},
  title = {{Implicit Surface Model Generator}},
  version = {1.0.0},
  year = {2021}
}
```

# Description

![Screenshot of the model generator...](/ss01.png)

The basic, single surface algorithm generates a volume in R3 and uses Marching Cubes to generate the mesh.
To offset the surfaces, I used the method described by Geier and merged them.



# References

Geier, M., & Alihussein, H. (2021). Computation of implicit representation of volumetric shells with predefined thickness. Algorithms, 14(4). https://doi.org/10.3390/a14040125
