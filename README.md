# nestcpp

An effort to implement 2D irregual bin packing for laser cutters and cnc milling. 2D IBP is about packing shapes (irregual polygons) as tightly as possbile to safe material when cutting from sheets. This project started out as a porting attempt of [Svgnest](http://svgnest.com) but after some initial work i decided to implement it from scratch.
A first step has been done by implementing no-fit polygon generation in a [C++ library](https://github.com/kallaballa/libnfp) which is a fundamental optimization needed for 2D IBP.
