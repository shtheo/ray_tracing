# Génération d'image par lancé de rayons

Théo Ponton - MOS 2.2 Infographie - Ecole Centrale de Lyon

## Introduction

Le but de ce projet est de générer des images à l'aide de lancé de rayons. Il s'agit de la même méthode qui est utilisé dans le cinéma afin de réaliser des images de synthèse. Ce code a été écrit en C++ afin d'avoir une compilation des plus rapides. Voici un exemple de résultat final : 

![100R_avant_transfo_geom](Release\100R_avant_transfo_geom.png)

Nous allons procéder étape par étape afin de montrer chaque fonctionnalité que permet le code. Voici les différentes étapes que nous traiterons :

- Création d'une image en C++ 
- Affichage d'une sphère seule
- Affichage de plusieurs sphères / notion de scène
- Surfaces spéculaire et transparentes
- Correction gamma
- Anti-Aliasing : adoucir les contours
- Eclairage indirect
- Lumière étendue / Ombres douces
- Affichage d'un maillage 3D 
- Les textures



## I - Algorithme général

L'image que nous allons créer est le résultat d'intersections multiples entre des rayons (objects Ray) lancés depuis une grille de pixels et une objet Scene. Cette Scene est composée d'un certain nombre d'objets Object :

- Des sphères (objet Sphere)
- Des triangles (objet Triangle)
- Des maillages constituées de plusieurs petits triangles (objet Geometry)

Nous avons donc une caméra, constituée d'un centre C et d'un angle d'ouverture fov.

```cpp
// Camera
double fov = 60 * M_PI / 180.;
Vector C(0, 0, 55);
```

Depuis cette caméra, nous allons lancer des rayons (Ray) à travers un grille de pixels. Un rayon est composé d'un point de départ attribut nommé C et d'un vecteur directeur u, attribut nommé u. Nous en lancerons un certain nombre par pixel. Si non lançons par exemple 30 rayons par pixel, nous allons lancer au total 512 * 512 * 30 = 7, 9 millions de rayons lancés.

```cpp
// Creating the pixels
for (int i = 0; i < W; i++) { // horizontally
    for (int j = 0; j < H; j++) { // vertically
        Vector pixelColor;
        for (int k=0; k<100; k++)
        {
            Vector CPrime = ...
            Vector uPrime = ...
            uPrime.normalize();
            Ray r = Ray(CPrime, uPrime);
            pixelColor = pixelColor + scene.getColor(r, epsilon, 5, n_sphere, n_air);
```

Comme nous pouvons le voir, pour chaque pixel (i,j), le but est d'avoir la couleur que dois avoir le pixel grâce à la méthode getColor() de la classe Scene. 

Pour cela, l'objet Scene contient une liste d'Object. Elle va donc dans la méthode getColor() utiliser sa méthode intersect afin de savoir si le rayon va toucher ou non un des ses objets et si oui lequel (le plus près du rayon). 

Une fois l'objet de l'intersection déterminée, il faut analyser sa surface. Si elle est transparente, spéculaire, diffuse ou si c'est la lumière. Selon la surface, le choix du rayon de lumière qui va continuer son chemin sera différent. La couleur du pixel en question sera la somme de émissivité de la surface touchée et de la fonction getColor() appelée sur le rayon réfléchi ou réfracté. Le nombre de rebond est majoré (ici il sera de 5 maximum pour ne pas trop perdre en performance).

```cpp
Vector getColor(Ray& r, double& epsilon, int bounce, double n_sphere, double n_air)
{
    Vector pixelColor(0, 0, 0);
    ...
    if (intersect(r, P, N, S, color) && bounce > 0) {
        if (S->mirror)
        {
            ...
            pixelColor = getColor(r_reflect, epsilon, bounce - 1, n_sphere, n_air);
        }
        else if (S->transparent)
        {
            double dotP = dot(r.u, N);
            if (dotP > 0) {
                ...
                pixelColor = getColor(r_ref, epsilon, bounce - 1, n_sphere, n_air);
            }
            else 
            {
                ...
                pixelColor = getColor(r_ref, epsilon, bounce - 1, n_sphere, n_air);
            }
        }
        else if (S->light)
        {
            pixelColor = Vector(1, 1, 1);
        }
        else
        {
			...
            if (intersect(rprime, Prime, Nprime, Sprime, colorprime)) { 
                Vector PrimeP = P - Prime;
                if (PrimeP.norm2() < dist_lum)
                {
                    pixelColor = Vector(0, 0, 0);
                }
                else
                {
                    pixelColor = I / (4 * M_PI*dist_lum) * (color / M_PI) * 
                        dot(PxPrime, N) * dot(-1. * PxPrime, Nprimee) / dot(Nprimee, OX);
                }
            }
            else
            {
                pixelColor = I / (4 * M_PI*dist_lum) * (color / M_PI) * dot(PxPrime, N) * 
                    dot(-1. * PxPrime, Nprimee) / dot(Nprimee, OX);
            }
            pixelColor = pixelColor + color * getColor(r_diffusion, epsilon, bounce, n_sphere, n_air);
        }
    }
    return pixelColor;
};
```

