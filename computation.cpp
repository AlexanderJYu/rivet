/**********************************************************************
Copyright 2014-2016 The RIVET Developers. See the COPYRIGHT file at
the top-level directory of this distribution.

This file is part of RIVET.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/

#include "computation.h"
#include "dcel/arrangement.h"
#include "dcel/arrangement_builder.h"
#include "debug.h"
#include "math/multi_betti.h"
#include "timer.h"
#include <chrono>
//#include "computing_s.h"
//#include <QtCore>
//#include <QApplication>
//#include <QString>


Computation::Computation(InputParameters& params, Progress& progress)
    : params(params)
    , progress(progress)
    , verbosity(params.verbosity)
{
}

Computation::~Computation()
{
}

std::unique_ptr<ComputationResult> Computation::compute_raw(ComputationInput& input)
{
    if (verbosity >= 2) {
        debug() << "\nBIFILTRATION:";
        debug() << "   Number of simplices of dimension " << params.dim << " : " << input.bifiltration_data().get_size(params.dim);
        debug() << "   Number of simplices of dimension " << (params.dim + 1) << " : " << input.bifiltration_data().get_size(params.dim + 1);
        if (verbosity >= 4) {
            debug() << "   Number of x-exact:" << input.x_exact.size();
            if (input.x_exact.size()) {
                std::ostringstream oss;
                oss << "; values" << input.x_exact.front() << " to " << input.x_exact.back();
                //debug() << QString::fromStdString(oss.str());
                debug() << oss.str().c_str();
            }
            debug() << "   Number of y-exact:" << input.y_exact.size();
            if (input.y_exact.size()) {
                std::ostringstream oss;
                oss << "; values" << input.y_exact.front() << " to " << input.y_exact.back();
                //debug() << QString::fromStdString(oss.str());
                debug() << oss.str().c_str();
            }
        }
    }
    //STAGE 3: COMPUTE MULTIGRADED BETTI NUMBERS

    std::unique_ptr<ComputationResult> result(new ComputationResult);
    //compute xi_0 and xi_1 at all bigrades
    if (verbosity >= 2) {
        debug() << "COMPUTING xi_0, xi_1, AND xi_2 FOR HOMOLOGY DIMENSION " << params.dim << ":";
    }
    FIRep rep(input.bifiltration(), verbosity);
    MultiBetti mb(rep, params.dim);
    Timer timer;
    mb.compute(result->homology_dimensions, progress);
    mb.compute_xi2(result->homology_dimensions);

    if (verbosity >= 2) {
        debug() << "  -- xi_i computation took " << timer.elapsed() << " milliseconds";
    }

    //store the xi support points
    mb.store_support_points(result->template_points);

    //store the dendrogram support points
    //SimplexTree* bif = &input.bifiltration();
    BifiltrationData* bif = &input.bifiltration_data();
    std::cout << "========================= bifiltration: " << std::endl;
    bif->print_bifiltration();
    computing_s* cs = new computing_s(bif,bif->num_x_grades(),bif->num_y_grades(),input.x_exact,input.y_exact);
    cs->compute_s();
    cs->print_T_0();
    cs->print_T_1();
    std::vector<TemplatePoint> S_as_template_pts;
    cs->store_S_points(S_as_template_pts);
    result->dendrogram_template_points = S_as_template_pts;
    std::cout << "========================= store_S_points: " << std::endl;
    for (auto& point : S_as_template_pts)
    {
        std::cout << "(" << point.x << "," << point.y << ")";
    }
    std::cout << std::endl;
    auto oracle= cs->get_oracle();


    template_points_ready(TemplatePointsMessage{ input.x_label, input.y_label, result->template_points, result->homology_dimensions, input.x_exact, input.y_exact }); //signal that xi support points are ready for visualization
    dendrogram_template_points_ready(TemplatePointsMessage{ input.x_label, input.y_label, result->dendrogram_template_points,
                                                            result->homology_dimensions, input.x_exact, input.y_exact });
    progress.advanceProgressStage(); //update progress box to stage 4

    //STAGES 4 and 5: BUILD THE LINE ARRANGEMENT AND COMPUTE BARCODE TEMPLATES

    //build the arrangement
    if (verbosity >= 2) {
        debug() << "CALCULATING ANCHORS AND BUILDING THE DCEL ARRANGEMENT";
    }

    timer.restart();
    ArrangementBuilder builder(verbosity);
    auto arrangement = builder.build_arrangement(mb, input.x_exact, input.y_exact, result->template_points, progress); ///TODO: update this -- does not need to store list of xi support points in xi_support
    // build dendrogram_arrangement
    ArrangementBuilder dendrogram_builder(verbosity);
    SimplexTree dummy_tree(0,0);
    auto dendrogram_arrangement = dendrogram_builder.build_arrangement(dummy_tree, *bif, input.x_exact, input.y_exact, S_as_template_pts, progress, oracle, cs);
    //NOTE: this also computes and stores barcode templates in the arrangement

    if (verbosity >= 2) {
        debug() << "   building the line arrangement and computing all barcode templates took"
                << timer.elapsed() << "milliseconds";
    }

    //send (a pointer to) the arrangement back to the VisualizationWindow
    debug() << "about to handle arrangement_ready signal";
    arrangement_ready(arrangement);
    // send (a pointer to) the dendrogram arrangement back to the VisualizationWindow
    debug() << "about to handle dendrogram_arrangement_ready signal";
    dendrogram_arrangement_ready(dendrogram_arrangement);
    //re-send xi support and other anchors
    if (verbosity >= 8) {
        debug() << "Sending" << result->template_points.size() << "anchors";
    }
    template_points_ready(TemplatePointsMessage{ input.x_label, input.y_label, result->template_points,
                                                 result->homology_dimensions, input.x_exact, input.y_exact });

    dendrogram_template_points_ready(TemplatePointsMessage{ input.x_label, input.y_label, result->dendrogram_template_points,
                                                            result->homology_dimensions, input.x_exact, input.y_exact });
    if (verbosity >= 10) {
        //TODO: Make this a separate flag, rather than using verbosity?
        arrangement->test_consistency();
    }
    result->arrangement = std::move(arrangement);
    result->dendrogram_arrangement = std::move(dendrogram_arrangement);
    return result;
}

std::unique_ptr<ComputationResult> Computation::dendrogram_compute_raw(ComputationInput& input)
{
    if (verbosity >= 2) {
        debug() << "\nBIFILTRATION:";
        debug() << "   Number of simplices of dimension " << params.dim << " : " << input.bifiltration_data().get_size(params.dim);
        debug() << "   Number of simplices of dimension " << (params.dim + 1) << " : " << input.bifiltration_data().get_size(params.dim + 1);
        if (verbosity >= 4) {
            debug() << "   Number of x-exact:" << input.x_exact.size();
            if (input.x_exact.size()) {
                std::ostringstream oss;
                oss << "; values" << input.x_exact.front() << " to " << input.x_exact.back();
                //debug() << QString::fromStdString(oss.str());
                debug() << oss.str().c_str();
            }
            debug() << "   Number of y-exact:" << input.y_exact.size();
            if (input.y_exact.size()) {
                std::ostringstream oss;
                oss << "; values" << input.y_exact.front() << " to " << input.y_exact.back();
                //debug() << QString::fromStdString(oss.str());
                debug() << oss.str().c_str();
            }
        }
    }
    //STAGE 3: COMPUTE MULTIGRADED BETTI NUMBERS

    std::unique_ptr<ComputationResult> result(new ComputationResult);
    /*
    //compute xi_0 and xi_1 at all bigrades
    if (verbosity >= 2) {
        debug() << "COMPUTING xi_0, xi_1, AND xi_2 FOR HOMOLOGY DIMENSION " << params.dim << ":";
    }
    MultiBetti mb(input.bifiltration_data(), params.dim);
    Timer timer;
    mb.compute(result->homology_dimensions, progress);
    mb.compute_xi2(result->homology_dimensions);

    if (verbosity >= 2) {
        debug() << "  -- xi_i computation took " << timer.elapsed() << " milliseconds";
    }

    //store the xi support points
    mb.store_support_points(result->template_points);*/

    //store the dendrogram support points
    //SimplexTree* bif = &input.bifiltration();
    BifiltrationData* bif = &input.bifiltration_data();
    std::cout << "========================= bifiltration: " << std::endl;
    bif->print_bifiltration();
    computing_s* cs = new computing_s(bif,bif->num_x_grades(),bif->num_y_grades(),input.x_exact,input.y_exact);
    cs->compute_s();
    cs->print_T_0();
    cs->print_T_1();
    std::vector<TemplatePoint> S_as_template_pts;
    cs->store_S_points(S_as_template_pts);
    result->dendrogram_template_points = S_as_template_pts;
    std::cout << "========================= store_S_points: " << std::endl;
    for (auto& point : S_as_template_pts)
    {
        std::cout << "(" << point.x << "," << point.y << ")";
    }
    std::cout << std::endl;
    auto oracle= cs->get_oracle();


    //template_points_ready(TemplatePointsMessage{ input.x_label, input.y_label, result->template_points, result->homology_dimensions, input.x_exact, input.y_exact }); //signal that xi support points are ready for visualization
    dendrogram_template_points_ready(TemplatePointsMessage{ input.x_label, input.y_label, result->dendrogram_template_points,
                                                            result->homology_dimensions, input.x_exact, input.y_exact });
    progress.advanceProgressStage(); //update progress box to stage 4

    //STAGES 4 and 5: BUILD THE LINE ARRANGEMENT AND COMPUTE BARCODE TEMPLATES

    //build the arrangement
    if (verbosity >= 2) {
        debug() << "CALCULATING ANCHORS AND BUILDING THE DCEL ARRANGEMENT";
    }

    //timer.restart();
    //ArrangementBuilder builder(verbosity);
    //auto arrangement = builder.build_arrangement(mb, input.x_exact, input.y_exact, result->template_points, progress); ///TODO: update this -- does not need to store list of xi support points in xi_support
    // build dendrogram_arrangement
    ArrangementBuilder dendrogram_builder(verbosity);
    SimplexTree dummy_tree(0,0);
    auto dendrogram_arrangement = dendrogram_builder.build_arrangement(dummy_tree, *bif, input.x_exact, input.y_exact, S_as_template_pts, progress, oracle, cs);
    //NOTE: this also computes and stores barcode templates in the arrangement
    /*
    if (verbosity >= 2) {
        debug() << "   building the line arrangement and computing all barcode templates took"
                << timer.elapsed() << "milliseconds";
    }*/

    //send (a pointer to) the arrangement back to the VisualizationWindow
    //debug() << "about to handle arrangement_ready signal";
    //arrangement_ready(arrangement);
    // send (a pointer to) the dendrogram arrangement back to the VisualizationWindow
    debug() << "about to handle dendrogram_arrangement_ready signal";
    dendrogram_arrangement_ready(dendrogram_arrangement);
    //re-send xi support and other anchors
    /*if (verbosity >= 8) {
        debug() << "Sending" << result->template_points.size() << "anchors";
    }
    template_points_ready(TemplatePointsMessage{ input.x_label, input.y_label, result->template_points,
                                                 result->homology_dimensions, input.x_exact, input.y_exact });*/

    dendrogram_template_points_ready(TemplatePointsMessage{ input.x_label, input.y_label, result->dendrogram_template_points,
                                                            result->homology_dimensions, input.x_exact, input.y_exact });
    /*if (verbosity >= 10) {
        //TODO: Make this a separate flag, rather than using verbosity?
        arrangement->test_consistency();
    }
    result->arrangement = std::move(arrangement);*/
    result->dendrogram_arrangement = std::move(dendrogram_arrangement);
    return result;
}

std::unique_ptr<ComputationResult> Computation::compute(InputData data)
{
    progress.advanceProgressStage(); //update progress box to stage 3

    auto input = ComputationInput(data);
    //print bifiltration statistics
    //debug() << "Computing from raw data";
    // return compute_raw(input);
    return dendrogram_compute_raw(input);
}
