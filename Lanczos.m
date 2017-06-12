%% Sean Reiter, June 8, 2017
clc
clear all
n = 100; %int size for random nxn Hermitian matrix 
A = zeros(n, n); %preallocate nxn zero matrix to increase computational efficiency
for i = 1:n %initiate for loop to generate upper triangular portion of random nxn Hermitian matrix
    for j = 1:i %initiate for loop to iterate through row i of random nxn Hermitian matrix
        A(i,j) = complex(rand, rand); %select a random complex no. w/ real and imag parts in the interval (0,1) as the ij-th element of A
    end
end
A = ctranspose(A) + A; %add A to its adjoint to make A an nxn Hermitian matrix

u_1 = zeros(n, 1); %preallocate nx1 zero vector to increase computational efficiency
for i = 1:n %initiate loop to generate arbitrary nx1 complex column vector
    u_1(i,1) = complex(rand, rand); %select a random complex no. w/ real and imag parts in the interval (0,1) as the i-th element of u_1
end
v_1 = u_1/(norm(u_1)); %norm arbitrary vector to have euclidean norm = 1, set to be 1st basis vector for our Krylov space

m = 20; %set number of iterations for Lanczos to run
w = A*v_1;  
a = ctranspose(w)*v_1;
u = w - a*v_1; %orthogonalize u against first basis vector  v_1
T = zeros(m,m); %preallocate space for mxm tridiagonal matrix T s.t. T = V*AV where V is the nxm matrix with columns as the basis vectors for our Krylov space
T(1,1) = a; %set first diagonal element of T = a_1 = w_1A = v_1*Av_1 
for i = 2:m %initiate loop of m interations to generate T
    b = norm(u); %calculate off diagonal elements of T, b, as the euclidean norm of u, the previously orthogonalized vector
    v = u/b; %norm u which was orthogonalized against all other vectors in basis in the previous iteration, and set as resulting vector v as the next basis vector
    w = A*v;
    a = ctranspose(w)*v; %calculate i-th diagonal element of T, a
    u = w - a*v - b*v_1; %orthogonalize w against basis vectors generated thus far, u will be normed in the next iteration and set as the next basis vector
    T(i,i) = a; %set i-th diagonal element of T = a
    T(i-1,i) = b; %set off diagonal elements equal to b
    T(i,i-1) = T(i-1,i); 
    v_1 = v; %save v_1 for next iteration as current v, to be used in the calculation of b
end
eig(A) %calculate eigenvalues of A
eig(T) %calculate eigenvalues of T
