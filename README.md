# Simulador de sensor acustico distribuido Phi OTDR
Este simulador fue desarrollado por Marcos/Sasmosaw Rojas Mardones como parte de su proyecto de titulación para la obtención de grado de ingeniero civil eléctronico con mención en telecomunicaciones y submención en fisica

El mismo consiste en un simulador de un sensor acoustico distribuido en fibra óptica basado en el principio de reflectometria óptica en el dominio temporal sensible a la fase, o simplemente, Phase Sensitive OTDR.

## Descripción

El contenido de cada carpeta se describe a continuación:
* Libraries : Contiene todas los archivos que definen el simulador, para hacer un nuevo proyectos deben de importarse dichos archivos.

* Examples : Contiene ejemplos de uso del simulador, dichos ejemplos pueden ser usados para estudiar el funcionamiento del simulador, o usar dichos codigos como plantilla par así construir nuevas simulaciones y/o experimentos, recomendado pues estan comentados en su totalidad explicando cada linea de codigo.

* CodesUsedInMemory : Contiene todos los archivos utilizados durante el desarrollo de la memoria, son utiles para estudiar el funcionamiento del simulador, ver varios ejemplos de usos y otros, sin embargo dichos archivos no tienen porque estar comentados, lo que puede dificultar el uso de los mismos

Todas las clases y funciones tienen incluidas de manera nativa documentación con el comando ```help```, por ejemplo, para saber más del uso de la clase ```classRayleigh```, en la terminal se puede introducir:
```
help classRayleigh
```
Para buscar más información sobre metodos especificos de la clase se pueder hacer con:
```
help classRayleigh.metodoEspecifico
```
## Modo de uso (Inicio rapido)

Para hacer uso de este simulador simplemente hay copiar y pegar los archivos dentro del directorio "libraries" dentro del proyecto de interes.

Se recomienda encarecidamente revisar los codigos de ejemplo al estudiar este procedimiento

El uso recomendado del simulador es el siguiente (Y el usado en los ejemplos):
- Inicializar todos los parametros de entrada para la creación de los modulos o clases
- Crear las instancias de las clases invocando los constructores respectivos, en este punto es buena idea generar algunas variables que no son generadas de manera nativa en la creación de las clases, por ejemplo:
```
%Frecuencia lenta
fml = Propagador.Vp/(2*fiber_1.fiberLength);
tml = 1/fml;

%Factor de conversión termica
FcT=1/(thermoOpticCoefficient+thermalExpansionCoefficient*Propagador.n);
%Factor de conversion tensorial
epsilon = -(1/2)*(Propagador.n^2)*((1-poissonRatio)*p12 - poissonRatio*p11);
FcE = 1/(Propagador.n*(epsilon + 1));
```
Donde ```Propagador``` es una instancia de la clase ```classPragation``` y ```fiber_1``` es una instancia de la calse ```classFiber```
- Calibración del factor de calibración, este procedimiento es explicado a mayor profundidad en la respectiva memoria, pero dicho de manera simple, se genera un factor de correción para punto/s de la fibra de interes por medio de la media geometrica
- Simulación, se genera un vector temporal, y las perturbaciones con sus caracteristicas espacio temporales, no existe manera nativa para hacer variaciones dinamicas, ni temporal ni espacialmente, por lo cual ambas deben ser realizadas por iteraciones de variaciones estaticas, para más detalles revisar los códigos

- Luego de este punto se tienen todos los resultados de interes para hacer analísis de la simulación, para ver casos de ejemplo especificos revise los archivos en Examples y Libraries