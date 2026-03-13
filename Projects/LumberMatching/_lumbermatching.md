# Lumber Matching Project

This project was done in May 2025 as an exploratory study in lumber matching. The idea is to use as little timber to design a spanning system as possible by cutting reciprocal segments from prismatic beams. This addresses the issue of offcuts associated with shaped timber beams, where material that is cut away is typically downcycled, by designing the offcut to also have structural use elsewhere in the system.

The algorithm works as follows:

1. draw a number of moment and shear demands from some roof or floor design (powered by Asap InternalForces component)
2. associate the combined load envelope (piecewise function) with the minimum height in timber that would analytically be required to support that load
3. slice the minimum heights at key points (using a convex hull) to create a simpler polygon representation of the irregular shape that could be cut with a bandsaw and satisfies the maximum height requirements from the more complex distribution at all points
4. given a set of available timber dimensions (uniform heights), find the optimal distribution of reciprocal segments
5. plot!

One of the main limitations with this project is that it is not possible (presently) to effectively and cheaply manufacture moment connections in timber. While this is tragic, maybe the project will come in handy some other time.
