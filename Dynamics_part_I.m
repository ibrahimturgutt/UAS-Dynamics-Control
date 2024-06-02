syms lmbd
A = [0.5-410*lmbd, 410; -0.45, -401.8*lmbd];
determinanat = det(A) %164738*lmbd^2-(2009*lmbd)/10+369/2
vector_determinant = [164738, -2009/10, 369/2];
roots_of_vector = roots(vector_determinant)