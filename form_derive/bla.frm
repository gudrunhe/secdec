#Define maxdepth "10"

S U, F, s, t, eps;
AutoDeclare S z;
AutoDeclare CF V;
CF Sum, Product;

L foo = Sum(z1, Product(z2, V(z3), Sum(z0, 3*z0)), z0);

print;
.sort

S e1, e2, e3, e4;
CF Derivative;

L bar = Derivative(foo, z0);

.sort


repeat;

id Derivative(Sum(e1?), e4?) = Derivative(e1, e4);
id Derivative(Sum(e1?, ?e2), e4?) = Derivative(e1, e4) + Derivative(Sum(?e2), e4);

id Derivative(Product(e1?), e4?) = Derivative(e1, e4);
id Derivative(Product(e1?, ?e2), e4?) = Product(Derivative(e1, e4), Product(?e2)) + 
					Product(e1, Derivative(Product(?e2), e4));


id Derivative(e1?, e1?) = 1;

#do e = {z0,z1,z2,z3,z4}
    id Derivative(`e', e2?) = 0;
#enddo


#do index = {1,...,`maxdepth'}
#do otherindex = 1,`index'
Argument;
#enddo

id Derivative(Sum(e1?), e4?) = Derivative(e1, e4);
id Derivative(Sum(e1?, ?e2), e4?) = Derivative(e1, e4) + Derivative(Sum(?e2), e4);

id Derivative(Product(e1?), e4?) = Derivative(e1, e4);
id Derivative(Product(e1?, ?e2), e4?) = Product(Derivative(e1, e4), Product(?e2)) + 
					Product(e1, Derivative(Product(?e2), e4));


id Derivative(e1?, e1?) = 1;

#do e = {z0,z1,z2,z3,z4}
    id Derivative(`e', e2?) = 0;
#enddo

#do otherindex = 1,`index'
EndArgument;
#enddo

#enddo


endrepeat;

print +s;

.end
