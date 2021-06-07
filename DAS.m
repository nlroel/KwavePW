classdef DAS
    %DAS 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        kgrid
        transducer
        offsets
        sensor_data
        scan_lines
        c0
    end
    
    methods
        function obj = DAS(kgrid,transducer,c0)
            obj.kgrid = kgrid;
            obj.transducer = transducer;
            delays = obj.transducer.beamforming_delays();
            obj.offsets = max(abs(delays))+delays;
            obj.c0 = c0;
        end
        
        function obj = load(obj, sensor_data)
            obj.sensor_data = sensor_data;
        end
        
        function scan_lines = get_scan_lines(obj,d)
            scan_lines = obj.das_beamforming(d);
        end
        
        function sl = das_beamforming(obj, d)
            sl = zeros(obj.transducer.number_elements ,obj.kgrid.Nt);
            scan_lines_ext = zeros(obj.kgrid.Nt, obj.transducer.number_elements + d);
            d_2 = round(d/2);
            scan_lines_ext(:, d_2:d_2+obj.transducer.number_elements-1) = obj.sensor_data';
            wd = ones(d,1);
            wd_norm = wd / sum(wd(:));
            delay_ext = inf(obj.kgrid.Nt, obj.transducer.number_elements + d);
            for j = 1:obj.transducer.number_elements
                delay = obj.getDelay(j);
                delay_ext(:, d_2:d_2+obj.transducer.number_elements-1) = delay;
                tmp_scan_line = interpMat2(scan_lines_ext(:, j:j+d-1),delay_ext(:, j:j+d-1),'linear',0);
                sl(j,:) =  tmp_scan_line * wd_norm;
            end
        end
        
        function delay = getDelay(obj, iy1)
            [YY,ZZ] = meshgrid(1:obj.transducer.number_elements,1:obj.kgrid.Nt);
            sy = YY * obj.kgrid.dy;
            y1 = iy1 * obj.kgrid.dy;
            sz = ZZ * obj.kgrid.dt * obj.c0 / 2.0;
            delay_theta = sz + obj.offsets(iy1)*obj.kgrid.dt*obj.c0;
            delay_d = sqrt((sy-y1).^2+sz.^2);
            delay = (delay_theta+delay_d)./(obj.kgrid.dt*obj.c0);
        end
        
    end
end

