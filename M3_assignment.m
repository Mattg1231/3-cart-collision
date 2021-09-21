% Matthew Gjanci
% ECE 202 Exercise M3
% 9/13/2021
% Elastic collision between three cars in 1D
%Discription: In this assignment we are looking for the final outcome of
%a 1d elastic collision that happens between three carts. For my
%calculations, my first collision starts with cart 2 and 3

clear

% ------GIVEN VARIABLES------

vi = [36 9 -45]; % An array of the Velocities for carts 1,2,&3 respectively. 
%Units in cm/s
m = [240 60 300]; % An array of the masses for carts 1,2,&3 respectively. 
%Units in grams

pi = sum(vi.*m); % Initial momentum of the system used for checking 
%equations. Units g*cm/s
KE = sum(vi.^2.*m/2); %Initial Kinetic energy(KE) of the system used 
%for checking equations. Units erg

M12 = m(1)+m(2); %Combining masses for carts 1&2 and 2&3 respectively to 
M23 = m(2)+m(3); %help shorten the equations. Units grams



% -------CALCULATIONS-------
%collision A (carts 2 and 3)

% equations used, to find the final velocities of the carts 2&3 
% respectively after the collision, units in cm/s
vA(2) = (2*m(3)*vi(3) + (m(2)-m(3))*vi(2))/M23; 
vA(3) = (2*m(2)*vi(2) + (m(3)-m(2))*vi(3))/M23;
vA(1) = vi(1) %cart 1's velocity stays the same as previous collision since
%it did not collide in this collision
%checking velocities after collision A through command window

% ------- checking ------

%checking answer through conservation of momentum 
dpA = pi - sum(vA.*m) % change in momentum (delta p1) should be 0 if 
%conserved unless point float error. Units g*cm/s

%checking answer through conservation of energy
dKEA = KE - sum(vA.^2.*m/2) % change in Kenetic energy (delta KE) should be
%0 if conserved unless point float error. Units erg

%------CHECKING FOR MORE COLLISIONS------

if vA(1)>vA(2) || vA(3)<vA(2) 
%If the left most cart's(cart 1) velocity is greater in the positive 
%direction than the middle cart(cart 2), then there will be a collision on 
%the left since cart one will either catch up to cart two or the other way 
%around. Similarly, if the right most cart's velocity is less than the
%middle cart's in the postive direction, then there will be a collision on 
%the right for the same reasoning as before.
    disp("There is another collision.")	
    
%If the first statement fails, then there cant be a collision.  
else
    disp("There is no collision.")
end



%collision B (carts 1 and 2)

% equations used, to find the final velocities of the carts 1&2 
% respectively after the collision, units in cm/s
vB(1) = (2*m(2)*vA(2) + (m(1)-m(2))*vA(1))/M12; 
vB(2) = (2*m(1)*vA(1) + (m(2)-m(1))*vA(2))/M12;
vB(3) = vA(3) %cart 3's velocity stays the same as previous collision since
%it did not collide in this collision
%checking velocities after collision B through command window

% ------- checking ------

%checking answer through conservation of momentum 
dpB = pi - sum(vB.*m) % change in momentum (delta p1) should be 0 if 
%conserved unless point float error. Units g*cm/s

%checking answer through conservation of energy
dKEB = KE - sum(vB.^2.*m/2) % change in Kenetic energy (delta KE) should be
%0 if conserved unless point float error. Units erg

%------CHECKING FOR MORE COLLISIONS------

if vB(1)>vB(2) || vB(3)<vB(2) 
%If the left most cart's(cart 1) velocity is greater in the positive 
%direction than the middle cart(cart 2), then there will be a collision on 
%the left since cart one will either catch up to cart two or the other way 
%around. Similarly, if the right most cart's velocity is less than the
%middle cart's in the postive direction, then there will be a collision on 
%the right for the same reasoning as before.
    disp("There is another collision.")	
%If the first statement fails, then there cant be a collision. 
else
    disp("There is no collision.")
end



%collision C (carts 2 and 3)

% equations used, to find the final velocities of the carts 2&3 
% respectively after the collision, units in cm/s
vC(2) = (2*m(3)*vB(3) + (m(2)-m(3))*vB(2))/M23; 
vC(3) = (2*m(2)*vB(2) + (m(3)-m(2))*vB(3))/M23; 
vC(1) = vB(1) %cart 1's velocity stays the same as previous collision since
%it did not collide in this collision
%checking velocities after collision C through command window

% ------- checking ------

%checking answer through conservation of momentum 
dpC = pi - sum(vC.*m) % change in momentum (delta p1) should be 0 if 
%conserved unless point float error. Units g*cm/s

%checking answer through conservation of energy
dKEC = KE - sum(vC.^2.*m/2) % change in Kenetic energy (delta KE) should be
%0 if conserved unless point float error. Units erg


%------CHECKING FOR MORE COLLISIONS------

if vC(1)>vC(2) || vC(3)<vC(2) 
%If the left most cart's(cart 1) velocity is greater in the positive 
%direction than the middle cart(cart 2), then there will be a collision on 
%the left since cart one will either catch up to cart two or the other way 
%around. Similarly, if the right most cart's velocity is less than the
%middle cart's in the postive direction, then there will be a collision on 
%the right for the same reasoning as before.
    disp("There is another collision.")

%If the first statement fails, then there cant be a collision.  
else
	disp("There is no collision.")
end



%collision D (carts 1 and 2)

% equations used, to find the final velocities of the carts 1&2 
% respectively after the collision, units in cm/s
vD(1) = (2*m(2)*vC(2) + (m(1)-m(2))*vC(1))/M12; 
vD(2) = (2*m(1)*vC(1) + (m(2)-m(1))*vC(2))/M12; 
vD(3) = vC(3) %cart 3's velocity stays the same as previous collision since
%it did not collide in this collision
%checking velocities after collision D through command window

% ------- checking ------

%checking answer through conservation of momentum 
dpD = pi - sum(vD.*m) % change in momentum (delta p1) should be 0 if 
%conserved unless point float error. Units g*cm/s

%checking answer through conservation of energy
dKED = KE - sum(vD.^2.*m/2) % change in Kenetic energy (delta KE) should be
%0 if conserved unless point float error. Units erg


%------CHECKING FOR MORE COLLISIONS------

if vD(1)>vD(2) || vD(3)<vD(2) 
%If the left most cart's(cart 1) velocity is greater in the positive 
%direction than the middle cart(cart 2), then there will be a collision on 
%the left since cart one will either catch up to cart two or the other way 
%around. Similarly, if the right most cart's velocity is less than the
%middle cart's in the postive direction, then there will be a collision on 
%the right for the same reasoning as before.
    disp("There is another collision.")	
%If the first statement fails, then there cant be a collision.  
else
    disp("There is no collision.")    
end



%collision E (carts 2 and 3)

% equations used, to find the final velocities of the carts 2&3 
% respectively after the collision, units in cm/s
vE(2) = (2*m(3)*vD(3) + (m(2)-m(3))*vD(2))/M23;
vE(3) = (2*m(2)*vD(2) + (m(3)-m(2))*vD(3))/M23; 
vE(1) = vD(1) %cart 1's velocity stays the same as previous collision since
%it did not collide in this collision
%checking velocities after collision E through command window

% ------- checking ------

%checking answer through conservation of momentum 
dpE = pi - sum(vE.*m) % change in momentum (delta p1) should be 0 if 
%conserved unless point float error. Units g*cm/s

dKEE = KE - sum(vE.^2.*m/2) % change in Kenetic energy (delta KE) should be
%0 if conserved unless point float error. Units ergor


%------CHECKING FOR MORE COLLISIONS------

if vE(1)>vE(2) || vE(3)<vE(2) 
%If the left most cart's(cart 1) velocity is greater in the positive 
%direction than the middle cart(cart 2), then there will be a collision on 
%the left since cart one will either catch up to cart two or the other way 
%around. Similarly, if the right most cart's velocity is less than the
%middle cart's in the postive direction, then there will be a collision on 
%the right for the same reasoning as before.
    disp("There is another collision.")	
%If the first statement fails, then there cant be a collision.  
else
    disp("There is no collision.")  
end

%I am satisified with my results. I think I had one point float error on 
%collision C since the kenetic energy came out to a very small number it 
%leads me to believe it was just a point float error. Other than that all 
%my equations seem to satisfy the assignment.
%There were a total of 5 collisions in this assigment.