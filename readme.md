# SPAR
The SPAR (Summed Product Analytic Representation) framework, is designed to serve as a general black-box interface through which kinetic energy operators (KEOs), potential energy functions (PEFs), dipole moment functions (DMFs) in sum-of-products form and possibly other operators can be communicated to a varational solver for ro-vibrational calculations of polyatomic molecules. 

SPAR is used to represent operators derived in sum-of-products form:

$`h(\mathbf{q}) = 
\sum_{{l}_1,{l}_2,\ldots,l_{3\mathcal{N}-6}}a_{{l}_1,{l}_2,\ldots,l_{3\mathcal{N}-6}} f^{(1)}_{l_1}(q_1)  f^{(2)}_{l_2}(q_2) \cdots f^{(3\mathcal{N}-6)}_{l_{3\mathcal{N}-6}}(q_{3\mathcal{N}-6})  \equiv 
\sum_{\mathbf{l}} a_{\mathbf{l}} \prod^{3\mathcal{N} - 6}_{k = 1} f^{(k)}_{l_k}(q_k)`$