function z = nonlinear_tophat_3D(U,d,K)

U1 = box_average_3D(U,d);
U2 = box_average_3D(U,K*d);

z = U.*U1./(U2.*U2);




