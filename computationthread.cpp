#include "computationthread.h"

#define CEREAL_SERIALIZE_FUNCTION_NAME cerealize

#include "dcel/mesh.h"
#include "interface/input_manager.h"
#include "interface/input_parameters.h"
#include "math/multi_betti.h"
#include "math/simplex_tree.h"
#include "math/xi_point.h"
#include "dcel/serialization.h"

#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>
#include "interface/console_interaction.h"

#include <QDebug>
#include <QTime>
#include <QProcess>
#include <QString>


ComputationThread::ComputationThread(InputParameters& params, QObject *parent) :
    QThread(parent),
    params(params),
    xi_support(),
    x_exact(),
    y_exact()
{ }

ComputationThread::~ComputationThread()
{ }

void ComputationThread::compute()
{
    start();
}

//this function does the work
void ComputationThread::run()
{
    QStringList args;

    args << QString::fromStdString(params.fileName)
    << QString::fromStdString(params.outputFile)
    << "-H" << QString::number(params.dim)
    << "-x" << QString::number(params.x_bins)
    << "-y" << QString::number(params.y_bins)
    << "-V" << QString::number(params.verbosity);
    auto console = RivetConsoleApp::start(args);

    if (!console->waitForStarted()) {
        qDebug() << "Error launching rivet_console:" << RivetConsoleApp::errorMessage(console->error()) ;
        return;
    }

    bool reading_xi = false;
    bool reading_mesh = false;
    std::stringstream ss;
    while(console->canReadLine() || console->waitForReadyRead()) {
        QString line = console->readLine();
        qDebug().noquote() << "console: " << line;
        if (reading_xi) {
            if (line.startsWith("END XI")) {
//                    cereal::JSONInputArchive archive(ss);
//                    cereal::BinaryInputArchive archive(ss);
                {
                    qDebug() << "Buffer is:";
                    qDebug() << QString::fromStdString(ss.str());
                    cereal::XMLInputArchive archive(ss);
                    XiSupportMessage message;
                    archive(message);
                    xi_support = message.xi_support;
                    x_exact = message.x_exact;
                    y_exact = message.y_exact;
                }
                reading_xi = false;
                emit xiSupportReady();
            } else {
                ss << line.toStdString();
            }
        } else if (reading_mesh) {
            if (line.startsWith("END ARRANGEMENT")) {
//                    cereal::JSONInputArchive archive(ss);
//                    cereal::BinaryInputArchive archive(ss);
                {
                    cereal::XMLInputArchive archive(ss);
                    archive(arrangement);
                }
                qDebug() << "Mesh received: " << arrangement->x_exact.size() << " x " << arrangement->y_exact.size();
                emit arrangementReady(&*arrangement);
                return;
            } else {
                ss << line.toStdString();
            }

        } else if (line.startsWith("PROGRESS ")) {
            auto progress = line.mid(QString("PROGRESS ").length()).trimmed();
            qDebug() << "***Progress is: " << progress;
            setCurrentProgress(progress.toInt());
        } else if (line.startsWith("STAGE")) {
            emit advanceProgressStage();
        } else if (line.startsWith("XI")) {
            reading_xi = true;
            ss.clear();
        } else if (line.startsWith("ARRANGEMENT")) {
            reading_mesh = true;
            ss.clear();
        }
    }

    qDebug() << "Mesh was not delivered";

}//end run()

