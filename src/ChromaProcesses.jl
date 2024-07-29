module ChromaProcesses

    using Chain
    using CSV
    using DataFrames
    using DataFramesMeta
    using CairoMakie

    export
        NoISTD,
        process_TIS_to_data_plots!

        
    struct TIS end
    struct TIC end
    struct MS end
    struct Raw end
    struct Log end
    struct Log10 end
    struct ISTD end
    struct NoISTD 
        min
        max
    end
    struct Heatmap end
    struct LinePlot end 
    struct SVG end
    struct PNG end
    struct PDF end

    # Text input functions
    function data_description_identifier(::TIS)
        "TIS"
    end

    function data_description_identifier(::TIC)
        "Chromatogram_RT-Abund"
    end

    function get_retention_time_name()
        "RetentionTimeMin"
    end

    function data_description_identifier(::MS)
        "Chromatogram_MS"
    end

    function xlab()
        "Retention Time (min)"
    end


    function ylab(::MS,::Union{Raw,Log,Log10})
        "m/z"
    end

    function ylab(::TIC,::Raw)
        "Abundance"
    end


    # Computational functions
    function get_identifier_based_file(data_typ,fold_pat)
        (root,_,_) = walkdir(fold_pat)
        identifier = data_description_identifier(data_typ)
        filename = filter(x->contains(x,identifier),root[3])[1]
        joinpath(fold_pat,filename)
    end

    function remove_istd(data,istd_min,istd_max)
        column_names_present = [i ∈ ["RetentionTimeMin","Abundance"] for i in names(data)] |>  all
        !(column_names_present) && error("RetentionTimeMin and Abundance are not present. Check using correct dataset.")
        @chain data begin
            @rtransform :Abundance = 
                :RetentionTimeMin > istd_min && :RetentionTimeMin < istd_max ? 
                0.0 : :Abundance
        end
    end

    function remove_istd(::TIC,istd_typ::NoISTD,data)
        column_names_present = [i ∈ ["RetentionTimeMin","Abundance"] for i in names(data)] |>  all
        !(column_names_present) && error("RetentionTimeMin and Abundance are not present. Check using correct dataset.")
        istd_min = istd_typ.min
        istd_max = istd_typ.max
        @chain data begin
            @rtransform :Abundance = 
                :RetentionTimeMin > istd_min && :RetentionTimeMin < istd_max ? 
                0.0 : :Abundance
        end
    end

    function remove_istd(::MS,istd_typ::NoISTD,data)
        istd_min = istd_typ.min
        istd_max = istd_typ.max
        istd_ind = findall(x-> (istd_min <= x <= istd_max),data.RetentionTimeMin)
        
        for col in Base.OneTo(size(data)[2]-1)
            col = col+1
            data[istd_ind,col] .= 0.0
        end
        data
    end

    function remove_istd(_,::ISTD,data)
        data
    end

    function transform_OpenChrom_data(::MS,path)
        @chain path begin
            CSV.read(_,DataFrame,header=true)
            DataFrames.rename(_, 2 => get_retention_time_name())
            select(_,Not([1,3]))
        end   
    end

    function transform_OpenChrom_data(::TIC,path)
        df = @chain path begin
            CSV.read(_,DataFrame,header=true)
            DataFrames.rename(_, 2 => get_retention_time_name())
        end

        DataFrame(
            RetentionTimeMin = df.RetentionTimeMin,
            Abundance = (sum ∘ eachcol)(df[:,4:end]))
    end

    function extract_data_from_TIS(::TIC,data_dir)
        path = get_identifier_based_file(TIS(),data_dir)
        df = transform_OpenChrom_data(TIC(),path)
        TIS_identifier = data_description_identifier(TIS())
        TIC_identifier = data_description_identifier(TIC())
        df_name = replace(basename(path), TIS_identifier => TIC_identifier)
        new_path = joinpath(data_dir,df_name)
        CSV.write(new_path,df,header=true)    
    end

    function extract_data_from_TIS(::MS,data_dir)
        path = get_identifier_based_file(TIS(),data_dir)
        df = transform_OpenChrom_data(MS(),path)
        TIS_identifier = data_description_identifier(TIS())
        MS_identifier = data_description_identifier(MS())
        df_name = replace(basename(path), TIS_identifier => MS_identifier)
        new_path = joinpath(data_dir,df_name)
        CSV.write(new_path,df,header=true)    
    end

    function get_plottable_data(::TIC,::Raw,df)
        (df.RetentionTimeMin,df.Abundance)
    end

    function get_plottable_data(::MS,::Raw,df)
        ions = parse.(Int,names(df)[2:end])
        rt = df[:,1]
        mdf = (Matrix(df[:,2:end]))
        (rt,ions,mdf)
    end

    function get_plottable_data(::MS,::Log,df)
        ions = parse.(Int,names(df)[2:end])
        rt = df[:,1]
        mdf = (Matrix(df[:,2:end]))
        lmdf = log.(mdf .+ 1)
        (rt,ions,lmdf)    
    end

    function get_plottable_data(::MS,::Log10,df)
        ions = parse.(Int,names(df)[2:end])
        rt = df[:,1]
        mdf = (Matrix(df[:,2:end]))
        lmdf = log10.(mdf .+ 1)
        (rt,ions,lmdf)    
    end

    function get_plottable_data(data_typ,tran_typ,istd_typ,data)
        @chain data begin
            remove_istd(data_typ,istd_typ,_)
            get_plottable_data(data_typ,tran_typ,_)
        end
    end

    function type_lowercase_string(T)
        @chain T begin
            string(_)
            replace(_,"()" => "")
            lowercase(_)
        end
    end

    function type_to_string(T)
        @chain T begin
            string(_)
            replace(_, r"[(].*$"=>"")
        end
    end

    function get_title_from_path(path)
        @chain path begin
            basename(_) 
            replace(_, "_" => " ") 
            replace(_, "-" => " ") 
            replace(_, ".csv" => "")
            replace(_, "060" => "6.0")
            replace(_, "105" => "10.5")
            replace(_, "pH" => "\npH")
        end
    end

    function plot_path_prefix(data_type,istd_typ,tran_typ,plot_type)
        plot_det = (data_type,istd_typ,tran_typ,plot_type)
        map(x->type_to_string(x),plot_det) |> x -> join(x,"_")
    end

    function plot_extension(plot_ext)
        type_lowercase_string(plot_ext)
    end

    function save_plot_path(data_typ,istd_typ,tran_typ,plot_typ,plot_ext,path)
        plot_pre = plot_path_prefix(data_typ,istd_typ,tran_typ,plot_typ)
        data_des = data_description_identifier(data_typ)
        plot_ext = plot_extension(plot_ext)
        path_ext = @chain path begin
            basename(_)
            split(_,".")
            _[end]
        end
        @chain path begin
                replace(_,data_des => plot_pre)
                replace(_,path_ext => plot_ext)
        end
    end

    function plot_base(::MS,fig,axis,legend,data)
        f = fig
        ax = axis
        heatmap!(ax,data...)
        legend
        f
    end

    function plot_base(::TIC,fig,axis,label,data)
        f = fig
        ax = axis
        lines!(ax,data...,label=label)
        axislegend()
        f
    end

    function generate_axis(placement,xlabel,ylabel,title)
        Axis(placement,xlabel = xlabel, ylabel = ylabel, title = title)
    end

    function get_axis(data_typ,tran_typ,path,placement)    
        title = get_title_from_path(path)
        ylabel = ylab(data_typ,tran_typ)
        generate_axis(placement,xlab(),ylabel,title)
    end

    function generate_fig()
        Figure(size=(1000,800),fontsize=30)
    end

    function generate_legend(::MS,tran_typ,data,placement)
        lab = string(type_to_string(tran_typ)," Intensity")
        Colorbar(placement,
            limits = 
            (minimum(data[end]),maximum(data[end])),
            vertical = false,
            label = lab)
    end

    function generate_legend(::TIC,tran_typ)
        string(type_to_string(tran_typ)," Abundance")
    end

    function plot_data(data_typ::MS,tran_typ,path,data)
        fig = generate_fig()
        ax = get_axis(data_typ,tran_typ,path,fig[1,1])    
        legend = generate_legend(data_typ,tran_typ,data,fig[2,1])
        plot_base(data_typ,fig,ax,legend,data)
    end

    function get_plot_type(::MS)
        Heatmap()
    end

    function get_plot_type(::TIC)
        LinePlot()
    end

    function plot_data(data_typ::TIC,tran_typ,path,data)
        fig = generate_fig()
        ax = get_axis(data_typ,tran_typ,path,fig[1,1])    
        label = generate_legend(data_typ,tran_typ)
        plot_base(data_typ,fig,ax,label,data)
    end

    function plot_save(data_typ,tran_typ,istd_typ,plot_ext,fold_pat)
        path = get_identifier_based_file(data_typ,fold_pat)
        df = CSV.read(path,DataFrame,header=true)
        data = get_plottable_data(data_typ,tran_typ,istd_typ,df)
        plot = plot_data(data_typ,tran_typ,path,data)
        plot_typ = get_plot_type(data_typ)
        save_path = save_plot_path(data_typ,istd_typ,tran_typ,plot_typ,plot_ext,path)
        save(save_path,plot)
    end

    function process_TIS_to_data_plots!(data_dir,no_istd)
        
        extract_data_from_TIS(TIC(),data_dir)
        extract_data_from_TIS(MS(),data_dir)
        plot_save(TIC(),Raw(),ISTD(),PNG(),data_dir)
        plot_save(TIC(),Raw(),ISTD(),PDF(),data_dir)
        plot_save(TIC(),Raw(),no_istd,PNG(),data_dir)
        plot_save(TIC(),Raw(),no_istd,PDF(),data_dir)
        plot_save(MS(),Raw(),ISTD(),PNG(),data_dir)
        plot_save(MS(),Raw(),ISTD(),PDF(),data_dir)
        plot_save(MS(),Raw(),no_istd,PNG(),data_dir)
        plot_save(MS(),Raw(),no_istd,PDF(),data_dir)
        plot_save(MS(),Log(),ISTD(),PNG(),data_dir)
        plot_save(MS(),Log(),ISTD(),PDF(),data_dir)
        plot_save(MS(),Log(),no_istd,PNG(),data_dir)
        plot_save(MS(),Log(),no_istd,PDF(),data_dir)
        plot_save(MS(),Log10(),ISTD(),PNG(),data_dir)
        plot_save(MS(),Log10(),ISTD(),PDF(),data_dir)
        plot_save(MS(),Log10(),no_istd,PNG(),data_dir)
        plot_save(MS(),Log10(),no_istd,PDF(),data_dir)
    end



end
