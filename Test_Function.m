function F = Test_Function( XX, D, testnumber )
%TEST FUNCTIONS
switch testnumber
    
    case 1                                                                 % Schwefel¡¯s problem 2.22 x [-10, 10]
        SUM1 = 0;
        SUM2 = 1;
        for i = 1:D
            SUM1 = SUM1 + abs(XX(i));
            SUM2 = SUM2 * abs(XX(i));
        end
        F = SUM1 + SUM2;
    
end
        
end