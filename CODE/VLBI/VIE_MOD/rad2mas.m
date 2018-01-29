function mas=rad2mas(rad)
% Converts radians to milliseconds of arc.
% Useage:  mas=rad2mas(rad)
% Input:   rad - angle in radians
% Output:  mas  - angle in milliseconds of arc

mas=rad.*180*3600/pi*1000;
