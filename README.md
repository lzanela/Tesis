## Repositorio para código de búsqueda de contraejemplos de Selection Monotonicity de distintos algoritmos de apportionment

- brewer contiene la implementación del algoritmo de Brewer, junto con código para buscar contraejemplos de selection monotonicity para dicho algoritmo.
- lexicographical_distribution contiene el problema lineal para buscar distribuciones que maximicen la probabilidad por orden lexicográfico, junto con código para buscar contraejemplos de selection monotonicity para dicho algoritmo.
- linear_optimal contiene el problema lineal para buscar distribuciones con vectores de costos genéricos, junto con código para buscar contraejemplos de selection monotonicity para dicho algoritmo.
- sampford contiene la implementación del algoritmo de Sampford, junto con código para buscar contraejemplos de strong selection monotonicity para dicho algoritmo.
- max_entropy contiene la implementación del algoritmo de máxima entropía, junto con código para buscar contraejemplos de selection monotonicity para dicho algoritmo.
- tesis_latex contiene el código de LaTex de la tesis.


#### Notas
- Probar contrarrecíproco para Sel. Mon. de Sampford 
- Agregar contraejemplo de sistema con quota + coallition incentive en Latex
- Hacer simplex a mano para n=4 y k=2 para los problemas con p y p', y ver si se puede establecer condiciones sobre p y eps para que las sols. no cumplan Sel. Mon.
- Bertsimas, Notes and sources 3.10: (Dantzig, 1963), Simplex para matrices sparse
- problema lineal con c = [39, 292, 712, 780, 147, 982, 473, 378, 687, 715] no encontró contraejemplo en un rato largo


#### Preguntas
- ¿La solución óptima de un problema lineal obtenida con el algoritmo Simplex depende del orden de las variables?
- Revisar preguntas del cuaderno
- Unequal probability sampling without replacement through a splitting method: ¿cuál es la motivación detrás de samplear los y_k con probas proporcionales a x_k?


#### Ideas
- Fractional Hedonic Games may serve as a way of killing all the coalition problems: given an allocation of the seats, one could design a utility function for each party over the other parties, based on the utility they would obtain by making a coalition with each of them. Those utility functions could define an unfeasible problem for defining coalitions.
- 