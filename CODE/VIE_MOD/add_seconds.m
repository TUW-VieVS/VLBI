function [yr, mm, dd, hh, minute, sec, doy] = add_seconds(yr, mm, dd, hh, minute, sec, dt_sec)
    % dt_sec kann positiv oder negativ sein
    % Schritt 1: alles in Sekunden seit Mitternacht des Tages umrechnen
    sec_total = hh*3600 + minute*60 + sec + dt_sec;

    % Schritt 2: Überträge auf Tage, Stunden, Minuten, Sekunden berechnen
    dd_shift = floor(sec_total / 86400);
    sec_total = sec_total - dd_shift*86400;
    if sec_total < 0
        sec_total = sec_total + 86400;
        dd_shift = dd_shift - 1;
    end

    hh = floor(sec_total / 3600);
    sec_total = sec_total - hh*3600;

    minute = floor(sec_total / 60);
    sec = sec_total - minute*60;

    % Schritt 3: Tag/Monat/Jahr anpassen
    dd = dd + dd_shift;
    while true
        % Tage im Monat berechnen (Schaltjahr berücksichtigen)
        if mm == 2
            if mod(yr,400)==0 || (mod(yr,4)==0 && mod(yr,100)~=0)
                days_in_month = 29;
            else
                days_in_month = 28;
            end
        elseif any(mm == [1,3,5,7,8,10,12])
            days_in_month = 31;
        else
            days_in_month = 30;
        end

        if dd > days_in_month
            dd = dd - days_in_month;
            mm = mm + 1;
            if mm > 12
                mm = 1;
                yr = yr + 1;
            end
        elseif dd < 1
            mm = mm - 1;
            if mm < 1
                mm = 12;
                yr = yr - 1;
            end
            % Tage des neuen Monats
            if mm == 2
                if mod(yr,400)==0 || (mod(yr,4)==0 && mod(yr,100)~=0)
                    days_in_month = 29;
                else
                    days_in_month = 28;
                end
            elseif any(mm == [1,3,5,7,8,10,12])
                days_in_month = 31;
            else
                days_in_month = 30;
            end
            dd = dd + days_in_month;
        else
            break
        end
    end

    dt = datetime([yr, mm, dd, hh, minute, sec]);
    doy = day(dt,'dayofyear');
end
