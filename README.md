
# mechanical linkage and FEM Softbody simulation in MATLAB
* fem2d.m   - softbody simulation using finite element method 
* linkages.m  - solving mechanical linkages differential equation

Open these files in MATLAB folder and type the following command on the command window 

## FEM 

* Neohookean
```
>fem2D(0)
```
![neo](/images/neo.gif)
* Saint Venant–Kirchhoff model
```
>fem2D(1)
```
![nvk](/images/nvk.gif)


## Linkage ()


* Intro
  
![0](/images/0.gif)

```
>linkages(0)
```

* drag-link
  
![1](/images/1.gif)
```
>linkages(1)
```

* double-rocker
  
![2](/images/2.gif)

```
>linkages(2)
```

* Hoekens
  
![3](/images/3.gif)

```
>linkages(3)
```

* Peaucellier-Lipkin
  
![4](/images/4.gif)

```
>linkages(4)
```

* Klann
  
![5](/images/5.gif)

```
>linkages(5)
```
