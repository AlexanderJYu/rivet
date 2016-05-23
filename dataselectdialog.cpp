#include "dataselectdialog.h"
#include "ui_dataselectdialog.h"

#include "interface/file_input_reader.h"
#include "interface/input_parameters.h"

#include <QDebug>
#include <QFileDialog>
#include <QFileInfo>
#include <fstream>


DataSelectDialog::DataSelectDialog(InputParameters& params, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DataSelectDialog),
    params(params),
    data_selected(false)
{
    ui->setupUi(this);

    ui->xbinSpinBox->setSpecialValueText(tr("No bins"));
    ui->ybinSpinBox->setSpecialValueText(tr("No bins"));

    //set initial values
    if(!params.fileName.empty())
        detect_file_type();
    ui->homDimSpinBox->setValue(params.dim);
    ui->xbinSpinBox->setValue(params.x_bins);
    ui->ybinSpinBox->setValue(params.y_bins);
}

DataSelectDialog::~DataSelectDialog()
{
    delete ui;
}

void DataSelectDialog::closeEvent(QCloseEvent* event)
{
    event->accept();

    if(!data_selected)
        qobject_cast<QWidget *>(this->parent())->close();
}

void DataSelectDialog::on_computeButton_clicked()
{
    params.dim = ui->homDimSpinBox->value();
    params.x_bins = ui->xbinSpinBox->value();
    params.y_bins = ui->ybinSpinBox->value();

    data_selected = true;

    emit dataSelected();

    close();
}

void DataSelectDialog::on_openFileButton_clicked()
{
    //prompt user to select a file
    QString selected_file = QFileDialog::getOpenFileName(this, tr("Open Data File"), "/ima/home/mlwright/Repos","All files (*.*);;Text files (*.txt)");

    if(!selected_file.isNull())
    {
      params.fileName = selected_file.toUtf8().constData();
      detect_file_type();
    }
}//end on_openFileButton_clicked()

void DataSelectDialog::detect_file_type()
{
  std::ifstream infile(params.fileName);

  if(!infile.is_open()) {
    invalid_file("Unable to read file.");
    return;
  }

  FileInputReader reader(infile);
  if(reader.has_next_line()) {
    invalid_file("Empty file.");
    return;
  }

  auto line = reader.next_line();
  auto file_type = line[0];

  auto supported = rivet::get_supported_type(file_type);

  if (!supported.first) {
    invalid_file("File type not recognized.");
    return;
  }

  auto supported_type = supported.second;
  ui->fileTypeLabel->setText("This file appears to contain " + QString::fromStdString(supported_type.description) + ".");
  QFileInfo fileInfo(QString::fromStdString(params.fileName));
  ui->fileLabel->setText("Selected file: " + fileInfo.fileName());

  //TODO: this updating of the params will need to happen in console also, need to refactor
  params.shortName = fileInfo.fileName().toUtf8().constData();
  params.raw_data = supported_type.is_data;

  ui->parameterFrame->setEnabled(params.raw_data);
  ui->computeButton->setEnabled(true);

}//end detect_file_type()

void DataSelectDialog::invalid_file(const QString &message) {
  ui->parameterFrame->setEnabled(false);
  ui->computeButton->setEnabled(false);
  ui->fileTypeLabel->setText(message);
}
